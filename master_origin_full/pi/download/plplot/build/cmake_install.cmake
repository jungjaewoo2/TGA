# Install script for directory: /home/pi/download/plplot/plplot-5.13.0

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
   "/home/pi/plplot/install_directory/share/doc/plplot/ABOUT;/home/pi/plplot/install_directory/share/doc/plplot/AUTHORS;/home/pi/plplot/install_directory/share/doc/plplot/COPYING.LIB;/home/pi/plplot/install_directory/share/doc/plplot/ChangeLog.release;/home/pi/plplot/install_directory/share/doc/plplot/Copyright;/home/pi/plplot/install_directory/share/doc/plplot/FAQ;/home/pi/plplot/install_directory/share/doc/plplot/NEWS;/home/pi/plplot/install_directory/share/doc/plplot/PROBLEMS;/home/pi/plplot/install_directory/share/doc/plplot/README;/home/pi/plplot/install_directory/share/doc/plplot/README.release;/home/pi/plplot/install_directory/share/doc/plplot/README.cumulated_release;/home/pi/plplot/install_directory/share/doc/plplot/README.testing")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/pi/plplot/install_directory/share/doc/plplot" TYPE FILE FILES
    "/home/pi/download/plplot/plplot-5.13.0/ABOUT"
    "/home/pi/download/plplot/plplot-5.13.0/AUTHORS"
    "/home/pi/download/plplot/plplot-5.13.0/COPYING.LIB"
    "/home/pi/download/plplot/plplot-5.13.0/ChangeLog.release"
    "/home/pi/download/plplot/plplot-5.13.0/Copyright"
    "/home/pi/download/plplot/plplot-5.13.0/FAQ"
    "/home/pi/download/plplot/plplot-5.13.0/NEWS"
    "/home/pi/download/plplot/plplot-5.13.0/PROBLEMS"
    "/home/pi/download/plplot/plplot-5.13.0/README"
    "/home/pi/download/plplot/plplot-5.13.0/README.release"
    "/home/pi/download/plplot/plplot-5.13.0/README.cumulated_release"
    "/home/pi/download/plplot/plplot-5.13.0/README.testing"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/pi/download/plplot/build/fonts/cmake_install.cmake")
  include("/home/pi/download/plplot/build/lib/cmake_install.cmake")
  include("/home/pi/download/plplot/build/include/cmake_install.cmake")
  include("/home/pi/download/plplot/build/src/cmake_install.cmake")
  include("/home/pi/download/plplot/build/data/cmake_install.cmake")
  include("/home/pi/download/plplot/build/bindings/cmake_install.cmake")
  include("/home/pi/download/plplot/build/drivers/cmake_install.cmake")
  include("/home/pi/download/plplot/build/utils/cmake_install.cmake")
  include("/home/pi/download/plplot/build/plplot_test/cmake_install.cmake")
  include("/home/pi/download/plplot/build/examples/cmake_install.cmake")
  include("/home/pi/download/plplot/build/scripts/cmake_install.cmake")
  include("/home/pi/download/plplot/build/doc/cmake_install.cmake")
  include("/home/pi/download/plplot/build/www/cmake_install.cmake")
  include("/home/pi/download/plplot/build/pkgcfg/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/pi/download/plplot/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
