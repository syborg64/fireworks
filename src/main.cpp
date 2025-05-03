#define GLM_FORCE_SWIZZLE
#define GLM_ENABLE_EXPERIMENTAL

#include <chrono>
#include <cstdlib>
#include <random>
#include <vector>

#include <glm/geometric.hpp>
#include <glm/fwd.hpp>

#include "AssetManager.hpp"
#include "Color.hpp"
#include "DeferredRenderer.hpp"
#include "Editor.hpp"
#include "Mesh.hpp"
#include "PbrMaterial.hpp"
#include "Physics.hpp"
#include "Resource.hpp"
#include "Stage.hpp"
#include "Time.hpp"
#include "collision/Collider.hpp"
#include "commands/EntityCommands.hpp"
#include "ecs.hpp"
#include "ecs/App.hpp"
#include "ecs/Entity.hpp"
#include "ecs/Query.hpp"
#include "ecs/commands/Commands.hpp"
#include "engine/Camera.hpp"
#include "engine/EnginePlugin.hpp"
#include "engine/Timer.hpp"
#include "engine/Transform.hpp"
#include "engine/Velocity.hpp"
#include "state.hpp"


using namespace cevy;
using namespace ecs;
using namespace engine;

const float base_mass = 20;

struct Focus {};

struct Particle {
  Time::time_point spawn_point;
  Time::time_point lifetime;
};

struct Fireburst {
  glm::vec3 color;
  float base_size;
  uint32_t recursions;
};

struct Sparkler {
  glm::vec3 color;
  float base_size;
  float half_life;
};

struct Trail {};

static glm::vec3 hsv2rgb(glm::vec3 c) {
  glm::vec4 K = glm::vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  glm::vec3 p = abs(fract(c.xxx() + K.xyz()) * 6.0f - K.www());
  return c.z * mix(K.xxx(), clamp(p - K.xxx(), 0.0f, 1.0f), c.y);
}

void initial_setup(Resource<asset::AssetManager> asset_manager, Commands cmd,
                   Resource<physics::Gravity> gravity, Resource<Atmosphere> atmosphere) {
  // gravity->acceleration *= 10;
  auto plane_mesh = asset_manager->add(primitives::plane(200, 4, 4), "plane.mesh");
  asset_manager->add(primitives::sphere(2, 6, 3), "sphere.mesh");
  auto mat = asset_manager->add(PbrMaterial(), "default.material");

  auto dark_mat = Handle<PbrMaterial>::Clone(mat);
  dark_mat->diffuse = {0, 0, 0};
  cmd.spawn(Transform(), plane_mesh, dark_mat);
  auto focus = cmd.spawn(Transform(), Focus{});

  auto _camera =
      cmd.spawn(Parent{focus.id()}, Camera(), Transform(glm::vec3(0, -100, 50), glm::quat({glm::radians(90.f), 0, 0}),
                                    glm::vec3(1)));

  // auto plane_handle = asset_manager->add(primitives::plane(256, 16, 16), "plane.mesh");
  // cmd.spawn(plane_handle, mat, Color(0.8, 0.8, 1), Transform());

  // glm::vec3 pos = glm::vec3(0, -10, 30);
  // glm::quat rot = glm::quatLookAt(-glm::normalize(pos), {0, 0, 1});
  // Transform tm = Transform(pos, rot, glm::vec3(.5, .5, .5));
  // SunLight light = {{1.3, 1.2, 0.9}, 30, 30};
  // cmd.spawn(tm, light);
}

void firework_emitter(Commands cmd, Query<const Camera, const Transform> camera,
                      Query<const Focus, const Transform> hit_plane,
                      Resource<Time> time, Resource<asset::AssetManager> asset_manager,
                      Resource<input::ButtonInput<input::MouseButton>> mouse_buttons,
                      cevy::ecs::Resource<cevy::input::cursorPosition> cursorPosition,
                      cevy::ecs::Resource<Window> window) {
  static cevy::physics::Collider collider(physics::Shape::Plane({0, 0, 0}, {0, 1, 0}));
  static int once = false;
  static std::optional<Entity> o_ent = {};
  if (camera.size() == 0 && !once) {
    once = true;
    return;
  }
  auto o_camera = camera.get_single();
  auto o_plane = hit_plane.get_single();

  if (mouse_buttons->is_just_pressed(input::MouseButton::Left) && o_camera && o_plane) {
    auto [camera, cam_trans] = o_camera.value();
    auto [_, plane_trans] = o_plane.value();
    auto window_size = window->windowSize();
    glm::vec2 screen_space = {(float(cursorPosition->pos.x) / window_size.x) * 2.f - 1.f,
                              ((float(cursorPosition->pos.y) / window_size.y) * -2.f + 1.f) /
                                  window_size.x * window_size.y};
    // screen_space = glm::inverse(camera.projection) * glm::vec4(screen_space, 1, 0);
    glm::mat4 cam_tm = glm::mat4(cam_trans.get_world());
    cam_tm /= cam_tm[3][3];
    auto origin = cam_tm * glm::vec4(0, 0, 0, 1);
    auto direction = glm::inverse(camera.projection) * glm::vec4(screen_space, 1, 1);
    direction = cam_tm * glm::vec4(direction.xyz(), 0);
    origin /= origin.w;
    auto ray = physics::Ray{origin, glm::normalize(direction.xyz())};

    auto collision = collider.raycast(plane_trans, ray);

    if (!collision.hit) {
      std::cout << "NO HIT" << std::endl;
      return;
    }
    std::cout << "HIT" << std::endl;

    auto distance = collision.location - glm::vec3(0, 0, 0);
    // v^2 = 2gh;
    auto h = distance.z * 1.1;
    auto g = 9.81;
    auto z_vel = std::sqrt(2 * h * g);
    auto time_to_h = z_vel / g / 1.5;

    glm::vec3 velocity = {distance.x / time_to_h / 2, distance.y / time_to_h / 2, z_vel};
    glm::vec3 spawn_point = {distance.x / 2, distance.y / 2, 0};
    Time::time_point lifetime = time->now();
    auto color = glm::vec3((std::rand() & 65535) / 65536.f, (std::rand() & 65535) / 65536.f,
                           (std::rand() & 65535) / 65536.f);
    color /= std::min({color.x, color.y, color.z});
    lifetime += Time::duration(time_to_h);
    // cmd.spawn(
    //     Particle{time->now(), lifetime}, cevy::physics::RigidBody(base_mass),
    //     // physics::Collider(0, physics::Shape::Sphere({}, 1)),
    //     Transform({collision.location.x / 2, 0, 0}, glm::quat({0, 0, 0}), glm::vec3(1)),
    //     TransformVelocity(velocity), asset_manager->get<Mesh>("sphere.mesh").value(),
    //     // Handle<Mesh>::Clone(asset_manager->get<Mesh>("sphere.mesh").value()),
    //     Handle<PbrMaterial>::Clone(asset_manager->get<PbrMaterial>("default.material").value()),
    //     Sparkler{color, 0.2, 4});

    // cmd.spawn(Particle{ time->now(), lifetime}, cevy::physics::RigidBody(base_mass),
    //           // physics::Collider(0, physics::Shape::Sphere({}, 1)),
    //           Transform({collision.location.x / 2, 0, 0}, glm::quat({0, 0, 0}), glm::vec3(1)),
    //           TransformVelocity(velocity), asset_manager->get<Mesh>("sphere.mesh").value(),
    //           // Handle<Mesh>::Clone(asset_manager->get<Mesh>("sphere.mesh").value()),
    //           Handle<PbrMaterial>::Clone(asset_manager->get<PbrMaterial>("default.material").value()),
    //           Fireburst{color, 1, 1});

    cmd.spawn(
        Particle{time->now(), lifetime}, cevy::physics::RigidBody(base_mass),
        // physics::Collider(0, physics::Shape::Sphere({}, 1)),
        Transform(spawn_point, glm::quat({0, 0, 0}), glm::vec3(1)),
        TransformVelocity(velocity), asset_manager->get<Mesh>("sphere.mesh").value(),
        // Handle<Mesh>::Clone(asset_manager->get<Mesh>("sphere.mesh").value()),
        Handle<PbrMaterial>::Clone(asset_manager->get<PbrMaterial>("default.material").value()),
        Fireburst{color, 1, 1}, Sparkler{color, 0.2, 2});
  }
}

void firework_exploder(Resource<Time> time,
                       Query<Entity, Transform, TransformVelocity, Particle, Fireburst,
                             Handle<Mesh>, Handle<PbrMaterial>>
                           firebursts,
                       Query<Entity, physics::RigidBody, Transform, TransformVelocity, Particle,
                             Sparkler, Handle<Mesh>, Handle<PbrMaterial>>
                           sparklers,
                       Commands cmd) {
  static std::random_device rd{};
  static std::mt19937 gen{rd()};
  static std::normal_distribution d{0.f, 1.f};

  for (auto [en, tm, vel, p, fw, mesh, mat] : firebursts) {
    bool explode = false;
    explode |= time->now() > p.lifetime;
    // explode |= std::abs(vel.position.z) < 0.1;
    explode |= glm::length(vel.position) < 0.1;
    explode |= tm.position.z < 0;

    if (explode) {
      cmd.despawn(en);
      const auto fragments = 100;
      if (fw.recursions > 0) {
        if (d(gen) > 0) {
          std::cout << "doubling" << std::endl;
          glm::vec3 color = {fw.color.y, fw.color.z, fw.color.x};
          color *= color;
          auto ext = d(gen) * 0.25 + 0.5;
          cmd.spawn(Particle{p.spawn_point, time->now() + Time::duration(ext)},
                    cevy::physics::RigidBody(base_mass),
                    // physics::Collider(0, physics::Shape::Sphere({}, 1)),
                    tm, TransformVelocity(glm::vec3(0, 0, 9.81 * ext)), mesh,
                    Handle<PbrMaterial>::Clone(mat), Fireburst{color, fw.base_size, fw.recursions});
        }

        for (uint i = 0; i < fragments; ++i) {
          float velocity = tm.position.z * (d(gen) + 100.f) / 20;
          glm::vec3 direction = glm::normalize(glm::vec3(d(gen), d(gen), d(gen)));
          direction *= velocity;
          auto ent =
              cmd.spawn(Particle{time->now(), time->now() + Time::duration(800 / velocity)},
                        cevy::physics::RigidBody(base_mass * fw.base_size / fragments),
                        Transform(tm.position), physics::Collider(0, physics::Shape::Sphere({}, 1)),
                        TransformVelocity(vel.position * 2.f + direction), mesh,
                        // Handle<Mesh>::Clone(mesh),
                        Handle<PbrMaterial>::Clone(mat),
                        Fireburst{fw.color, fw.base_size / std::pow(float(fragments), 1.f / 3.f),
                                  fw.recursions - 1});
          if ((gen() - gen.min()) / double(gen.max() - gen.min()) < 3. / fragments) {
            ent.insert(Sparkler{fw.color, fw.base_size / std::pow(float(fragments), 1.f / 3.f), 0.5});
            ent.insert(TransformVelocity(direction + glm::vec3{0, 0, 9.81}));
          }
        }
      }
    }
  }
  for (auto [en, body, tm, vel, p, fw, mesh, mat] : sparklers) {
    if (time->now() > p.lifetime || tm.position.z < 0) {
      cmd.despawn(en);
      continue;
    }
    auto probability = 1 - std::exp(-fw.half_life * time->delta_seconds());

    auto r = double(gen() - gen.min()) / double(gen.max() - gen.min());
    if (fw.half_life > 0 && r > probability) {
      // std::cout << "SPARKLING" << std::endl;
      float velocity = (d(gen) + 3.f) / 2;
      // float velocity = tm.position.z * (d(gen) + 100.f) / 80;
      glm::vec3 direction = glm::normalize(glm::vec3(d(gen), d(gen), d(gen)));
      direction *= velocity;
      auto lifetime = (p.lifetime - p.spawn_point) * 0.95 + time->now();

      cmd.spawn(Particle{p.spawn_point, lifetime}, cevy::physics::RigidBody(body.mass() * 0.05),
                Transform(tm.position, glm::quat({0, 0, 0}), {0, 0, 0}),
                physics::Collider(0, physics::Shape::Sphere({}, 3)),
                TransformVelocity(vel.position / 4.f + direction), mesh,
                // Handle<Mesh>::Clone(mesh),
                Handle<PbrMaterial>::Clone(mat), Sparkler{fw.color, fw.base_size * 0.5f, 0});
      // body.iMass /= 0.95;
    }
  }
}

void firework_renderer(
    Resource<Time> time,
    Query<Transform, option<TransformVelocity>, Particle, Fireburst, Handle<PbrMaterial>> fireworks,
    Query<Transform, option<TransformVelocity>, Particle, Sparkler, Handle<PbrMaterial>>
        sparklers) {
  for (auto [tm, o_vel, p, fw, mat] : fireworks) {
    float remaining = (p.lifetime - time->now()) / (p.lifetime - p.spawn_point);
    float energy = o_vel ? std::min(glm::length(o_vel->position) / 10, 1.f) * remaining : remaining;
    tm.scale = glm::vec3(fw.base_size * energy);
    mat->emit = fw.color * energy;
  }
  for (auto [tm, o_vel, p, fw, mat] : sparklers) {
    float remaining = (p.lifetime - time->now()) / (p.lifetime - p.spawn_point);
    // float energy = o_vel ? std::min(glm::length(o_vel->position) / 10, 1.f) * remaining :
    // remaining;
    tm.scale = glm::vec3(fw.base_size * remaining);
    if (fw.half_life != 0) {
      tm.scale = glm::vec3(0);
    }
    mat->emit = fw.color * remaining;
  }
}

void firework_trailer(
    Resource<Time> time, Commands cmd,
    Query<Transform, Particle, Fireburst, Handle<PbrMaterial>, Handle<Mesh>> fireworks,
    Query<Entity, Particle, Trail> trails) {
  for (auto [tm, p, fw, mat, mesh] : fireworks) {
    if ((random() & 665535) * time->delta_seconds() > 2048.f / (fw.recursions)) {
      auto lifetime = (p.lifetime - time->now()) * 0.5 + time->now();
      std::cout << "TRAILED" << std::endl;
      auto transform = tm;
      transform.scale = glm::vec3(0.2);
      // cmd.spawn(Particle{ time->now(), lifetime},
      //               transform,
      //               mesh,
      //               mat,
      //               Trail{});
    }
  }

  std::vector<Entity> dead;

  for (auto [en, p, _t] : trails) {
    if (time->now() > p.lifetime)
      cmd.despawn(en);
  }
  for (auto &d : dead) {
  }
}

void catchall(Resource<Time> time, Commands cmd,
              Query<Entity, Particle, option<Transform>> particles) {
  for (auto [en, p, o_tm] : particles) {
    if (time->now() > p.lifetime + Time::duration(0.5) || (o_tm && o_tm->position.z < 0)) {
      cmd.despawn(en);
    }
  }
}

void move_camera(Resource<input::ButtonInput<input::KeyCode>> keyboard,
                 Query<Focus, Transform> cam_q, Resource<ecs::Time> time) {
  glm::vec3 direction = {0, 0, 0};
  float speed = 10;

  for (auto [_, transform] : cam_q) {
    if (keyboard->is_pressed(input::KeyCode::A)) {
      direction.x -= 1;
    }
    if (keyboard->is_pressed(input::KeyCode::D)) {
      direction.x += 1;
    }
    if (keyboard->is_pressed(input::KeyCode::Q)) {
      direction.y -= 1;
    }
    if (keyboard->is_pressed(input::KeyCode::E)) {
      direction.y += 1;
    }
    if (keyboard->is_pressed(input::KeyCode::W)) {
      direction.z -= 1;
    }
    if (keyboard->is_pressed(input::KeyCode::S)) {
      direction.z += 1;
    }
    float delta_time = time->raw().count();

    if (glm::length(direction) != 0) {
      transform.translateXYZ(transform.rotation * glm::normalize(direction) * speed * delta_time);
    }
  }
}

void rotate_camera(
                   Query<Focus, Transform> table_q,
                   Query<Camera, Transform> cam_q,
                   Resource<input::ButtonInput<input::MouseButton>> mouse_buttons,
                   cevy::ecs::EventReader<input::mouseMotion> mouse_motion_reader) {
  static glm::vec2 rotation = {0 * glm::pi<float>(), glm::pi<float>() * 0.3f};

  if (mouse_buttons->is_pressed(input::MouseButton::Right)) {
    // std::cout << "mousing!" << std::endl;
    for (const auto &mouse_motion : mouse_motion_reader) {
      if (mouse_motion.delta.has_value()) {
        rotation.x -= mouse_motion.delta.value().x * 0.005;
        rotation.y -= mouse_motion.delta.value().y * 0.005;
        rotation.y = glm::clamp(rotation.y, 0.f, glm::pi<float>());
        // std::cout << "motion!" << mouse_motion.delta.value().x << " " <<
        // (mouse_motion.delta.value().y) << std::endl;
      } else {
        throw std::runtime_error("mouse motion has no value!!");
      }
      auto xQuat = glm::quat({0., 0., rotation.x});
      auto yQuat = glm::quat({rotation.y, 0., 0.});
      for (auto [_, transform] : cam_q) {
        transform.rotation = yQuat;
      }
      for (auto [_, transform] : table_q) {
        transform.rotation = xQuat;
      }
    }
  }
}

int main() {
  App app;
  app.init_component<Focus>();
  app.init_component<Particle>();
  app.init_component<Fireburst>();
  app.init_component<Sparkler>();
  app.init_component<Trail>();
  app.add_plugins(
      Engine<glWindow::Builder<cevy::engine::DeferredRenderer /* , editor::Editor */>>());
  app.add_plugins(physics::PhysicsPlugin());
  app.add_systems<core_stage::Startup>(initial_setup);
  app.add_systems<core_stage::Update>(rotate_camera);
  app.add_systems<core_stage::Update>(move_camera);
  app.add_systems<core_stage::Update>(firework_emitter);
  app.add_systems<core_stage::Update>(firework_exploder);
  app.add_systems<core_stage::Update>(catchall);
  app.add_systems<core_stage::PostUpdate>(firework_renderer);
  // app.add_systems<core_stage::Update>(firework_trailer);

  app.run();
  return 0;
}
