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
include drivers/CMakeFiles/mem.dir/depend.make

# Include the progress variables for this target.
include drivers/CMakeFiles/mem.dir/progress.make

# Include the compile flags for this target's objects.
include drivers/CMakeFiles/mem.dir/flags.make

drivers/CMakeFiles/mem.dir/mem.c.o: drivers/CMakeFiles/mem.dir/flags.make
drivers/CMakeFiles/mem.dir/mem.c.o: /home/pi/download/plplot/plplot-5.13.0/drivers/mem.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/download/plplot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object drivers/CMakeFiles/mem.dir/mem.c.o"
	cd /home/pi/download/plplot/build/drivers && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/mem.dir/mem.c.o   -c /home/pi/download/plplot/plplot-5.13.0/drivers/mem.c

drivers/CMakeFiles/mem.dir/mem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mem.dir/mem.c.i"
	cd /home/pi/download/plplot/build/drivers && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pi/download/plplot/plplot-5.13.0/drivers/mem.c > CMakeFiles/mem.dir/mem.c.i

drivers/CMakeFiles/mem.dir/mem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mem.dir/mem.c.s"
	cd /home/pi/download/plplot/build/drivers && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pi/download/plplot/plplot-5.13.0/drivers/mem.c -o CMakeFiles/mem.dir/mem.c.s

drivers/CMakeFiles/mem.dir/mem.c.o.requires:

.PHONY : drivers/CMakeFiles/mem.dir/mem.c.o.requires

drivers/CMakeFiles/mem.dir/mem.c.o.provides: drivers/CMakeFiles/mem.dir/mem.c.o.requires
	$(MAKE) -f drivers/CMakeFiles/mem.dir/build.make drivers/CMakeFiles/mem.dir/mem.c.o.provides.build
.PHONY : drivers/CMakeFiles/mem.dir/mem.c.o.provides

drivers/CMakeFiles/mem.dir/mem.c.o.provides.build: drivers/CMakeFiles/mem.dir/mem.c.o


# Object files for target mem
mem_OBJECTS = \
"CMakeFiles/mem.dir/mem.c.o"

# External object files for target mem
mem_EXTERNAL_OBJECTS =

drivers/mem.so: drivers/CMakeFiles/mem.dir/mem.c.o
drivers/mem.so: drivers/CMakeFiles/mem.dir/build.make
drivers/mem.so: src/libplplot.so.15.0.0
drivers/mem.so: /usr/lib/arm-linux-gnueabihf/libm.so
drivers/mem.so: drivers/CMakeFiles/mem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/download/plplot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C shared module mem.so"
	cd /home/pi/download/plplot/build/drivers && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mem.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
drivers/CMakeFiles/mem.dir/build: drivers/mem.so

.PHONY : drivers/CMakeFiles/mem.dir/build

drivers/CMakeFiles/mem.dir/requires: drivers/CMakeFiles/mem.dir/mem.c.o.requires

.PHONY : drivers/CMakeFiles/mem.dir/requires

drivers/CMakeFiles/mem.dir/clean:
	cd /home/pi/download/plplot/build/drivers && $(CMAKE_COMMAND) -P CMakeFiles/mem.dir/cmake_clean.cmake
.PHONY : drivers/CMakeFiles/mem.dir/clean

drivers/CMakeFiles/mem.dir/depend:
	cd /home/pi/download/plplot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/download/plplot/plplot-5.13.0 /home/pi/download/plplot/plplot-5.13.0/drivers /home/pi/download/plplot/build /home/pi/download/plplot/build/drivers /home/pi/download/plplot/build/drivers/CMakeFiles/mem.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : drivers/CMakeFiles/mem.dir/depend

