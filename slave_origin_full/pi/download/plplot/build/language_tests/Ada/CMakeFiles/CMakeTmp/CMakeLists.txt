
cmake_minimum_required(VERSION 3.6.2 FATAL_ERROR)
set(CMAKE_MODULE_PATH "/home/pi/download/plplot/build/language_tests/Ada/cmake/modules/language_support/cmake")
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "" FORCE)
project(test Ada)
add_executable(testadacompiler testadacompiler.adb)
  