# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/include

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
   "/home/pi/plplot/install_directory/include/plplot/disptab.h;/home/pi/plplot/install_directory/include/plplot/drivers.h;/home/pi/plplot/install_directory/include/plplot/pdf.h;/home/pi/plplot/install_directory/include/plplot/pldebug.h;/home/pi/plplot/install_directory/include/plplot/pldll.h;/home/pi/plplot/install_directory/include/plplot/plevent.h;/home/pi/plplot/install_directory/include/plplot/plplot.h;/home/pi/plplot/install_directory/include/plplot/plplotP.h;/home/pi/plplot/install_directory/include/plplot/plstrm.h;/home/pi/plplot/install_directory/include/plplot/plxwd.h;/home/pi/plplot/install_directory/include/plplot/plConfig.h;/home/pi/plplot/install_directory/include/plplot/plDevs.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/include/plplot" TYPE FILE FILES
    "/home/pi/download/plplot/plplot-5.13.0/include/disptab.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/drivers.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/pdf.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/pldebug.h"
    "/home/pi/download/plplot/build/include/pldll.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/plevent.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/plplot.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/plplotP.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/plstrm.h"
    "/home/pi/download/plplot/plplot-5.13.0/include/plxwd.h"
    "/home/pi/download/plplot/build/include/plConfig.h"
    "/home/pi/download/plplot/build/include/plDevs.h"
    )
endif()

