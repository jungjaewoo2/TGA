# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/examples/c++

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
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/plc++demos.h;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x01cc.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x00.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x01.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x02.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x03.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x04.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x05.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x06.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x07.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x08.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x09.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x10.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x11.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x12.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x13.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x14.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x15.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x16.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x17.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x18.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x19.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x20.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x21.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x22.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x23.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x24.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x25.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x26.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x27.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x28.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x29.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x30.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x31.cc;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/x33.cc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++" TYPE FILE FILES
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/plc++demos.h"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x01cc.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x00.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x01.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x02.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x03.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x04.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x05.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x06.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x07.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x08.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x09.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x10.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x11.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x12.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x13.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x14.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x15.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x16.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x17.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x18.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x19.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x20.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x21.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x22.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x23.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x24.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x25.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x26.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x27.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x28.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x29.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x30.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x31.cc"
    "/home/pi/download/plplot/plplot-5.13.0/examples/c++/x33.cc"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++" TYPE FILE FILES "/home/pi/download/plplot/plplot-5.13.0/examples/c++/CMakeLists.txt")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/c++" TYPE FILE RENAME "Makefile" FILES "/home/pi/download/plplot/build/examples/c++/Makefile.examples")
endif()

