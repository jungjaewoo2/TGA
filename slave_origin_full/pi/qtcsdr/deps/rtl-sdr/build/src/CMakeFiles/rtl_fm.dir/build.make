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
CMAKE_SOURCE_DIR = /home/pi/qtcsdr/deps/rtl-sdr

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pi/qtcsdr/deps/rtl-sdr/build

# Include any dependencies generated for this target.
include src/CMakeFiles/rtl_fm.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/rtl_fm.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/rtl_fm.dir/flags.make

src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o: src/CMakeFiles/rtl_fm.dir/flags.make
src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o: ../src/rtl_fm.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pi/qtcsdr/deps/rtl-sdr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o"
	cd /home/pi/qtcsdr/deps/rtl-sdr/build/src && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/rtl_fm.dir/rtl_fm.c.o   -c /home/pi/qtcsdr/deps/rtl-sdr/src/rtl_fm.c

src/CMakeFiles/rtl_fm.dir/rtl_fm.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/rtl_fm.dir/rtl_fm.c.i"
	cd /home/pi/qtcsdr/deps/rtl-sdr/build/src && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pi/qtcsdr/deps/rtl-sdr/src/rtl_fm.c > CMakeFiles/rtl_fm.dir/rtl_fm.c.i

src/CMakeFiles/rtl_fm.dir/rtl_fm.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/rtl_fm.dir/rtl_fm.c.s"
	cd /home/pi/qtcsdr/deps/rtl-sdr/build/src && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pi/qtcsdr/deps/rtl-sdr/src/rtl_fm.c -o CMakeFiles/rtl_fm.dir/rtl_fm.c.s

src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.requires:

.PHONY : src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.requires

src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.provides: src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.requires
	$(MAKE) -f src/CMakeFiles/rtl_fm.dir/build.make src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.provides.build
.PHONY : src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.provides

src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.provides.build: src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o


# Object files for target rtl_fm
rtl_fm_OBJECTS = \
"CMakeFiles/rtl_fm.dir/rtl_fm.c.o"

# External object files for target rtl_fm
rtl_fm_EXTERNAL_OBJECTS =

src/rtl_fm: src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o
src/rtl_fm: src/CMakeFiles/rtl_fm.dir/build.make
src/rtl_fm: src/libconvenience_static.a
src/rtl_fm: src/librtlsdr.so.0.5git
src/rtl_fm: /usr/lib/arm-linux-gnueabihf/libusb-1.0.so
src/rtl_fm: src/CMakeFiles/rtl_fm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pi/qtcsdr/deps/rtl-sdr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable rtl_fm"
	cd /home/pi/qtcsdr/deps/rtl-sdr/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rtl_fm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/rtl_fm.dir/build: src/rtl_fm

.PHONY : src/CMakeFiles/rtl_fm.dir/build

src/CMakeFiles/rtl_fm.dir/requires: src/CMakeFiles/rtl_fm.dir/rtl_fm.c.o.requires

.PHONY : src/CMakeFiles/rtl_fm.dir/requires

src/CMakeFiles/rtl_fm.dir/clean:
	cd /home/pi/qtcsdr/deps/rtl-sdr/build/src && $(CMAKE_COMMAND) -P CMakeFiles/rtl_fm.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/rtl_fm.dir/clean

src/CMakeFiles/rtl_fm.dir/depend:
	cd /home/pi/qtcsdr/deps/rtl-sdr/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pi/qtcsdr/deps/rtl-sdr /home/pi/qtcsdr/deps/rtl-sdr/src /home/pi/qtcsdr/deps/rtl-sdr/build /home/pi/qtcsdr/deps/rtl-sdr/build/src /home/pi/qtcsdr/deps/rtl-sdr/build/src/CMakeFiles/rtl_fm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/rtl_fm.dir/depend

