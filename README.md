# Download
to download this project, first download the source code ZIP from the `<> Code` button

Then, on the same page, click on the `Cevy @ xxxx` submodule in the folders in files to go to the Cevy project's ghithub

Then, download the Cevy source ZIP using the same `<> Code` button. Make sure you didn't change the version before downloading. The string after the name should match between the download and what's referenced in the Fireworks project

Once you've downloaded and extract both ZIPs, simply move the Cevy ZIP contents into the empty Cevy folder in Fireworks;

Your folder tree should end up something like this:
```
Fireworks/
├─ Cevy/
│  ├─ CMakeLists.txt
│  ├─ examples/
│  ├─ assets/
│  ├─ src/
├─ CMakeLists.txt
├─ README.md
├─ src/
```

# Clone
use the following command to clone with the Cevy submodule:
```
git clone --recurse-submodules git@github.com:Syborg64/Fireworks.git
```
if you've already cloned, you can run `git submodule init` and `git submodule update` to retrieve the submodule

note that the submodule uses ssh for cloning.

If you can't use ssh:
- clone the repo normally
- edit the `.gitmodules` file with the transport you wish to use
- run `git submodule init` and `git submodule update` to retrieve Cevy with your preferred transport

# Dependencies
This project has external and internal dependencies.

Please make sure the following libraries are installed and accessible to CMake:

for Windows:
- gl3w
- glfw3
- glm
- OpenGl

for linux:
- glew
- glfw3
- glm
- OpenGl

Cevy is bundled with internal dependencies, that are downloaded and built automatically:
- imgui

# Build
Firework is tested on windows to work with CMake, using ninja as a generator and VS 17 2022 as a compiler

on linux, it is tested with Unix Makefiles and clang

# Run
Make sure you run the project from the directory containing the binary. It needs the working directory to be set correctly to access engine assets like shader code

# Usage
Click on the screen to create a firework that will rise and explode where you clicked

Right-click lets you rotate the camera
