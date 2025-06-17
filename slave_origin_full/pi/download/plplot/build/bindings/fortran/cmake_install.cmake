# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/bindings/fortran

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
  foreach(file
      "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/libplplotfortran.so.0.0.0"
      "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/libplplotfortran.so.0"
      "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/libplplotfortran.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/home/pi/plplot/install_directory/lib")
    endif()
  endforeach()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/libplplotfortran.so.0.0.0;/home/pi/plplot/install_directory/lib/libplplotfortran.so.0;/home/pi/plplot/install_directory/lib/libplplotfortran.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib" TYPE SHARED_LIBRARY FILES
    "/home/pi/download/plplot/build/bindings/fortran/libplplotfortran.so.0.0.0"
    "/home/pi/download/plplot/build/bindings/fortran/libplplotfortran.so.0"
    "/home/pi/download/plplot/build/bindings/fortran/libplplotfortran.so"
    )
  foreach(file
      "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/libplplotfortran.so.0.0.0"
      "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/libplplotfortran.so.0"
      "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/libplplotfortran.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/pi/download/plplot/build/src:::"
           NEW_RPATH "/home/pi/plplot/install_directory/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plplot.mod;/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plplot_double.mod;/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plplot_single.mod;/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plplot_private_utilities.mod;/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plplot_graphics.mod;/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plplot_private_exposed.mod;/home/pi/plplot/install_directory/lib/fortran/modules/plplot/plplot_types.mod")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/fortran/modules/plplot" TYPE FILE FILES
    "/home/pi/download/plplot/build/bindings/fortran/plplot.mod"
    "/home/pi/download/plplot/build/bindings/fortran/plplot_double.mod"
    "/home/pi/download/plplot/build/bindings/fortran/plplot_single.mod"
    "/home/pi/download/plplot/build/bindings/fortran/plplot_private_utilities.mod"
    "/home/pi/download/plplot/build/bindings/fortran/plplot_graphics.mod"
    "/home/pi/download/plplot/build/bindings/fortran/plplot_private_exposed.mod"
    "/home/pi/download/plplot/build/bindings/fortran/plplot_types.mod"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/pkgconfig/plplot-fortran.pc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/pkgconfig" TYPE FILE FILES "/home/pi/download/plplot/build/pkgcfg/plplot-fortran.pc")
endif()

