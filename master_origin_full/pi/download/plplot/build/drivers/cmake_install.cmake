# Install script for directory: /home/pi/download/plplot/plplot-5.13.0/drivers

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
   "/home/pi/plplot/install_directory/share/doc/plplot/README.drivers;/home/pi/plplot/install_directory/share/doc/plplot/README.wxwidgets")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/doc/plplot" TYPE FILE FILES
    "/home/pi/download/plplot/plplot-5.13.0/drivers/README.drivers"
    "/home/pi/download/plplot/plplot-5.13.0/drivers/README.wxwidgets"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so"
         RPATH "/home/pi/plplot/install_directory/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE MODULE FILES "/home/pi/download/plplot/build/drivers/mem.so")
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so"
         OLD_RPATH "/home/pi/download/plplot/build/src:::"
         NEW_RPATH "/home/pi/plplot/install_directory/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.so")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/mem.driver_info")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE FILE FILES "/home/pi/download/plplot/build/drivers/mem.driver_info")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so"
         RPATH "/home/pi/plplot/install_directory/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE MODULE FILES "/home/pi/download/plplot/build/drivers/null.so")
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so"
         OLD_RPATH "/home/pi/download/plplot/build/src:::"
         NEW_RPATH "/home/pi/plplot/install_directory/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.so")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/null.driver_info")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE FILE FILES "/home/pi/download/plplot/build/drivers/null.driver_info")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so"
         RPATH "/home/pi/plplot/install_directory/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE MODULE FILES "/home/pi/download/plplot/build/drivers/ps.so")
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so"
         OLD_RPATH "/home/pi/download/plplot/build/src:::"
         NEW_RPATH "/home/pi/plplot/install_directory/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.so")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/ps.driver_info")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE FILE FILES "/home/pi/download/plplot/build/drivers/ps.driver_info")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so"
         RPATH "/home/pi/plplot/install_directory/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE MODULE FILES "/home/pi/download/plplot/build/drivers/svg.so")
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so"
         OLD_RPATH "/home/pi/download/plplot/build/src:::"
         NEW_RPATH "/home/pi/plplot/install_directory/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.so")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/svg.driver_info")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE FILE FILES "/home/pi/download/plplot/build/drivers/svg.driver_info")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so"
         RPATH "/home/pi/plplot/install_directory/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE MODULE FILES "/home/pi/download/plplot/build/drivers/xfig.so")
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so"
         OLD_RPATH "/home/pi/download/plplot/build/src:::"
         NEW_RPATH "/home/pi/plplot/install_directory/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.so")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xfig.driver_info")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE FILE FILES "/home/pi/download/plplot/build/drivers/xfig.driver_info")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so"
         RPATH "/home/pi/plplot/install_directory/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE MODULE FILES "/home/pi/download/plplot/build/drivers/xwin.so")
  if(EXISTS "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so"
         OLD_RPATH "/home/pi/download/plplot/build/src:::"
         NEW_RPATH "/home/pi/plplot/install_directory/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.so")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers/xwin.driver_info")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/lib/plplot5.13.0/drivers" TYPE FILE FILES "/home/pi/download/plplot/build/drivers/xwin.driver_info")
endif()

