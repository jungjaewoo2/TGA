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

# Include any dependencies generated for this target.
include drivers/CMakeFiles/svg.dir/depend.make

# Include the progress variables for this target.
include drivers/CMakeFiles/svg.dir/progress.make

# Include the compile flags for this target's objects.
include drivers/CMakeFiles/svg.dir/flags.make

drivers/CMakeFiles/svg.dir/svg.c.o: drivers/CMakeFiles/svg.dir/flags.make
drivers/CMakeFiles/svg.dir/svg.c.o: /home/pi/download/plplot/plplot-5.13.0/drivers/svg.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/download/plplot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object drivers/CMakeFiles/svg.dir/svg.c.o"
	cd /home/pi/download/plplot/build/drivers && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/svg.dir/svg.c.o   -c /home/pi/download/plplot/plplot-5.13.0/drivers/svg.c

drivers/CMakeFiles/svg.dir/svg.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/svg.dir/svg.c.i"
	cd /home/pi/download/plplot/build/drivers && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pi/download/plplot/plplot-5.13.0/drivers/svg.c > CMakeFiles/svg.dir/svg.c.i

drivers/CMakeFiles/svg.dir/svg.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/svg.dir/svg.c.s"
	cd /home/pi/download/plplot/build/drivers && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pi/download/plplot/plplot-5.13.0/drivers/svg.c -o CMakeFiles/svg.dir/svg.c.s

drivers/CMakeFiles/svg.dir/svg.c.o.requires:

.PHONY : drivers/CMakeFiles/svg.dir/svg.c.o.requires

drivers/CMakeFiles/svg.dir/svg.c.o.provides: drivers/CMakeFiles/svg.dir/svg.c.o.requires
	$(MAKE) -f drivers/CMakeFiles/svg.dir/build.make drivers/CMakeFiles/svg.dir/svg.c.o.provides.build
.PHONY : drivers/CMakeFiles/svg.dir/svg.c.o.provides

drivers/CMakeFiles/svg.dir/svg.c.o.provides.build: drivers/CMakeFiles/svg.dir/svg.c.o


# Object files for target svg
svg_OBJECTS = \
"CMakeFiles/svg.dir/svg.c.o"

# External object files for target svg
svg_EXTERNAL_OBJECTS =

drivers/svg.so: drivers/CMakeFiles/svg.dir/svg.c.o
drivers/svg.so: drivers/CMakeFiles/svg.dir/build.make
drivers/svg.so: src/libplplot.so.15.0.0
drivers/svg.so: /usr/lib/arm-linux-gnueabihf/libm.so
drivers/svg.so: drivers/CMakeFiles/svg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/download/plplot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C shared module svg.so"
	cd /home/pi/download/plplot/build/drivers && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/svg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
drivers/CMakeFiles/svg.dir/build: drivers/svg.so

.PHONY : drivers/CMakeFiles/svg.dir/build

drivers/CMakeFiles/svg.dir/requires: drivers/CMakeFiles/svg.dir/svg.c.o.requires

.PHONY : drivers/CMakeFiles/svg.dir/requires

drivers/CMakeFiles/svg.dir/clean:
	cd /home/pi/download/plplot/build/drivers && $(CMAKE_COMMAND) -P CMakeFiles/svg.dir/cmake_clean.cmake
.PHONY : drivers/CMakeFiles/svg.dir/clean

drivers/CMakeFiles/svg.dir/depend:
	cd /home/pi/download/plplot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/download/plplot/plplot-5.13.0 /home/pi/download/plplot/plplot-5.13.0/drivers /home/pi/download/plplot/build /home/pi/download/plplot/build/drivers /home/pi/download/plplot/build/drivers/CMakeFiles/svg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : drivers/CMakeFiles/svg.dir/depend

