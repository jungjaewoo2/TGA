# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/examples/c

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/pi/plplot/install_directory")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/plcdemos.h;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/tutor.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/test_plend.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/test_plbuf.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x00c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x01c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x02c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x03c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x04c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x05c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x06c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x07c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x08c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x09c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x10c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x11c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x12c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x13c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x14c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x15c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x16c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x17c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x18c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x19c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x20c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x21c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x22c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x23c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x24c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x25c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x26c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x27c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x28c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x29c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x30c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x31c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x32c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x33c.c;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/x34c.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c" TYPE FILE FILES
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/plcdemos.h"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/tutor.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/test_plend.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/test_plbuf.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x00c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x01c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x02c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x03c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x04c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x05c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x06c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x07c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x08c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x09c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x10c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x11c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x12c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x13c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x14c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x15c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x16c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x17c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x18c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x19c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x20c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x21c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x22c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x23c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x24c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x25c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x26c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x27c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x28c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x29c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x30c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x31c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x32c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x33c.c"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c/x34c.c"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c" TYPE FILE FILES "/home/pi/download/plplot/plplot-5.13.0/examples/c/CMakeLists.txt")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c" TYPE FILE RENAME "Makefile" FILES "/home/pi/download/plplot/build/examples/c/Makefile.examples")
endif()

