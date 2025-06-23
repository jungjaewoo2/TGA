cmake_minimum_required(VERSION 3.6.2 FATAL_ERROR)
project(test NONE)
# Locations where PLplot cmake build system first looks for cmake modules.
set(CMAKE_MODULE_PATH
  "/home/pi/download/plplot/plplot-5.13.0/cmake/modules"
  "/home/pi/download/plplot/plplot-5.13.0/cmake/modules/language_support/cmake"
  )

message(STATUS "CMAKE_GENERATOR = Unix Makefiles")
#enable_language(D OPTIONAL)
enable_language(D)
