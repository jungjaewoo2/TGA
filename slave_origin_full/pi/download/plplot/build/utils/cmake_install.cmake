# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/utils

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
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/bin/pltek" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/bin/pltek")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/bin/pltek"
         RPATH "/home/pi/plplot/install_directory/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/bin/pltek")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/bin" TYPE EXECUTABLE FILES "/home/pi/download/plplot/build/utils/pltek")
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/bin/pltek" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/bin/pltek")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/bin/pltek"
         OLD_RPATH "/home/pi/download/plplot/build/src:::"
         NEW_RPATH "/home/pi/plplot/install_directory/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/pi/plplot/install_directory/bin/pltek")
    endif()
  endif()
endif()

