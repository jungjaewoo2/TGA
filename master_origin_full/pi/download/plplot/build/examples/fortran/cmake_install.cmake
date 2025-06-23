# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/examples/fortran

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
   "/home/pi/plplot/install_directory/lib/libplfortrandemolib.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib" TYPE STATIC_LIBRARY FILES "/home/pi/download/plplot/build/examples/fortran/libplfortrandemolib.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plfortrandemolib.mod")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/fortran/modules/plplot" TYPE FILE FILES "/home/pi/download/plplot/build/examples/fortran/plfortrandemolib.mod")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x00f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x01f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x02f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x03f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x04f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x05f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x06f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x07f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x08f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x09f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x10f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x11f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x12f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x13f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x14f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x15f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x16f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x16af.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x17f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x18f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x19f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x20f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x21f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x22f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x23f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x24f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x25f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x26f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x27f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x28f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x29f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x30f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x31f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/x33f.f90;/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/README_precision")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran" TYPE FILE FILES
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x00f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x01f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x02f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x03f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x04f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x05f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x06f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x07f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x08f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x09f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x10f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x11f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x12f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x13f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x14f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x15f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x16f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x16af.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x17f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x18f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x19f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x20f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x21f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x22f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x23f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x24f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x25f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x26f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x27f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x28f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x29f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x30f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x31f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/x33f.f90"
    "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/README_precision"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran" TYPE FILE RENAME "Makefile" FILES "/home/pi/download/plplot/build/examples/fortran/Makefile.examples")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/plfortrandemos.inc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran" TYPE FILE FILES "/home/pi/download/plplot/build/examples/fortran/plfortrandemos.inc")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0/examples/fortran" TYPE FILE FILES "/home/pi/download/plplot/plplot-5.13.0/examples/fortran/CMakeLists.txt")
endif()

