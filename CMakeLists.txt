cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)

project(Cevy LANGUAGES CXX VERSION 0.1)
set(Cevy_VERSION 0.1)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (WIN32)
    set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ -lwsock32 -lws2_32 ${CMAKE_CXX_STANDARD_LIBRARIES}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive -Wl,-subsystem,windows -Wl,-allow-multiple-definition")
    #add_compile_options(-mwindows)
endif(WIN32)

add_subdirectory(Cevy)

# set(CMAKE_CXX_CLANG_TIDY clang-tidy -header-filter=. -checks=*,-llvmlibc-*)

add_subdirectory(src)

option(TESTS "Build tests" OFF)
if(TESTS)
    add_custom_target(tests-run
                        COMMENT "Executing tests"
                        COMMAND test-fireworks
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    )

    add_subdirectory(test)
endif()

option(RUN "Launch game" OFF)
if(RUN)
    add_custom_target(run
                        COMMENT "Launching binary"
                        COMMAND fireworks
                        # WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
                    )
    add_custom_target(run_debug
                        COMMENT "Launching binary with valgrind"
                        DEPENDS fireworks
                        # COMMAND pwd
                        COMMAND valgrind ./fireworks
                        # WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
                    )
endif()

option(SOUCRE_PACKAGE "Create source packaging" OFF)
if(SOURCE_PACKAGE)
    include("./cpack_source.cmake")
endif()
