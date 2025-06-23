# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/bindings

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/pi/download/plplot/build/bindings/swig-support/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/c++/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/fortran/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/tcl/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/tk/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/tk-x-plat/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/python/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/octave/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/java/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/wxwidgets/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/ada/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/d/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/ocaml/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/lua/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/qt_gui/cmake_install.cmake")

endif()

