# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/data

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
   "/home/pi/plplot/install_directory/share/plplot5.13.0/plstnd5.fnt;/home/pi/plplot/install_directory/share/plplot5.13.0/plxtnd5.fnt;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap0_default.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap0_alternate.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap0_white_bg.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap0_black_on_white.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap1_default.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap1_blue_yellow.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap1_blue_red.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap1_gray.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap1_radar.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap1_highfreq.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cmap1_lowfreq.pal;/home/pi/plplot/install_directory/share/plplot5.13.0/cglobe.map;/home/pi/plplot/install_directory/share/plplot5.13.0/globe.map;/home/pi/plplot/install_directory/share/plplot5.13.0/usa.map;/home/pi/plplot/install_directory/share/plplot5.13.0/usaglobe.map")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/plplot5.13.0" TYPE FILE FILES
    "/home/pi/download/plplot/plplot-5.13.0/data/plstnd5.fnt"
    "/home/pi/download/plplot/plplot-5.13.0/data/plxtnd5.fnt"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap0_default.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap0_alternate.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap0_white_bg.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap0_black_on_white.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap1_default.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap1_blue_yellow.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap1_blue_red.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap1_gray.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap1_radar.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap1_highfreq.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cmap1_lowfreq.pal"
    "/home/pi/download/plplot/plplot-5.13.0/data/cglobe.map"
    "/home/pi/download/plplot/plplot-5.13.0/data/globe.map"
    "/home/pi/download/plplot/plplot-5.13.0/data/usa.map"
    "/home/pi/download/plplot/plplot-5.13.0/data/usaglobe.map"
    )
endif()

