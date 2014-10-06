project(visualization.starburst)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

find_package(Kodi REQUIRED)
find_package(OpenGL REQUIRED)

add_definitions(-DHAS_SDL_OPENGL)
include_directories(${OpenGL_INCLUDE_DIR}
                    ${KODI_INCLUDE_DIR})

set(STARBURST_SOURCES src/StarBurst.cpp)

SET(DEPLIBS ${OPENGL_LIBRARIES})

build_addon(visualization.starburst STARBURST DEPLIBS)

include(CPack)
