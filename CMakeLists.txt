project(visualization.fountain)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

find_package(xbmc REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SOIL REQUIRED)

add_definitions(-DHAS_SDL_OPENGL)
include_directories(${OpenGL_INCLUDE_DIR}
                    ${GLEW_INCLUDE_DIR}
                    ${SOIL_INCLUDE_DIRS}
                    ${XBMC_INCLUDE_DIR})

set(FOUNTAIN_SOURCES src/Fountain.cpp
                     src/ParticleSystem.cpp
                     src/Util.cpp)

SET(DEPLIBS ${OPENGL_LIBRARIES}
             ${SOIL_LIBRARIES})

build_addon(visualization.fountain FOUNTAIN DEPLIBS)

include(CPack)
