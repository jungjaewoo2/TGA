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

# Utility rule file for test_dyndrivers.

# Include the progress variables for this target.
include drivers/CMakeFiles/test_dyndrivers.dir/progress.make

test_dyndrivers: drivers/CMakeFiles/test_dyndrivers.dir/build.make

.PHONY : test_dyndrivers

# Rule to build all files generated by this target.
drivers/CMakeFiles/test_dyndrivers.dir/build: test_dyndrivers

.PHONY : drivers/CMakeFiles/test_dyndrivers.dir/build

drivers/CMakeFiles/test_dyndrivers.dir/clean:
	cd /home/pi/download/plplot/build/drivers && $(CMAKE_COMMAND) -P CMakeFiles/test_dyndrivers.dir/cmake_clean.cmake
.PHONY : drivers/CMakeFiles/test_dyndrivers.dir/clean

drivers/CMakeFiles/test_dyndrivers.dir/depend:
	cd /home/pi/download/plplot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/download/plplot/plplot-5.13.0 /home/pi/download/plplot/plplot-5.13.0/drivers /home/pi/download/plplot/build /home/pi/download/plplot/build/drivers /home/pi/download/plplot/build/drivers/CMakeFiles/test_dyndrivers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : drivers/CMakeFiles/test_dyndrivers.dir/depend

