# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/pi/download/plplot/plplot-5.13.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pi/download/plplot/build

# Utility rule file for check_ifort_definitions_file.

# Include the progress variables for this target.
include bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/progress.make

bindings/fortran/CMakeFiles/check_ifort_definitions_file:
	cd /home/pi/download/plplot/build/bindings/fortran && /usr/bin/cmake -E echo Check\ that\ bindings/fortran/plplotfortran_ifort.def\ is\ consistent\ with\ the\ symbols\ in\ the\ plplotfortran\ library
	cd /home/pi/download/plplot/build/bindings/fortran && /usr/bin/cmake -E remove -f /home/pi/download/plplot/build/bindings/fortran/plplotfortran_ifort.def_compare
	cd /home/pi/download/plplot/build/bindings/fortran && bash /home/pi/download/plplot/build/bindings/fortran/plplotfortran_ifort.def.sh
	cd /home/pi/download/plplot/build/bindings/fortran && cmp /home/pi/download/plplot/plplot-5.13.0/bindings/fortran/plplotfortran_ifort.def /home/pi/download/plplot/build/bindings/fortran/plplotfortran_ifort.def_compare

check_ifort_definitions_file: bindings/fortran/CMakeFiles/check_ifort_definitions_file
check_ifort_definitions_file: bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/build.make

.PHONY : check_ifort_definitions_file

# Rule to build all files generated by this target.
bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/build: check_ifort_definitions_file

.PHONY : bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/build

bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/clean:
	cd /home/pi/download/plplot/build/bindings/fortran && $(CMAKE_COMMAND) -P CMakeFiles/check_ifort_definitions_file.dir/cmake_clean.cmake
.PHONY : bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/clean

bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/depend:
	cd /home/pi/download/plplot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/download/plplot/plplot-5.13.0 /home/pi/download/plplot/plplot-5.13.0/bindings/fortran /home/pi/download/plplot/build /home/pi/download/plplot/build/bindings/fortran /home/pi/download/plplot/build/bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : bindings/fortran/CMakeFiles/check_ifort_definitions_file.dir/depend

