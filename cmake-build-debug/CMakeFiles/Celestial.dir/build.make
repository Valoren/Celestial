# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /snap/clion/107/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/107/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/miquel/Desktop/Celestial

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/miquel/Desktop/Celestial/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Celestial.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Celestial.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Celestial.dir/flags.make

CMakeFiles/Celestial.dir/src/main.cpp.o: CMakeFiles/Celestial.dir/flags.make
CMakeFiles/Celestial.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/miquel/Desktop/Celestial/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Celestial.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Celestial.dir/src/main.cpp.o -c /home/miquel/Desktop/Celestial/src/main.cpp

CMakeFiles/Celestial.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Celestial.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/miquel/Desktop/Celestial/src/main.cpp > CMakeFiles/Celestial.dir/src/main.cpp.i

CMakeFiles/Celestial.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Celestial.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/miquel/Desktop/Celestial/src/main.cpp -o CMakeFiles/Celestial.dir/src/main.cpp.s

CMakeFiles/Celestial.dir/src/parser.cpp.o: CMakeFiles/Celestial.dir/flags.make
CMakeFiles/Celestial.dir/src/parser.cpp.o: ../src/parser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/miquel/Desktop/Celestial/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Celestial.dir/src/parser.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Celestial.dir/src/parser.cpp.o -c /home/miquel/Desktop/Celestial/src/parser.cpp

CMakeFiles/Celestial.dir/src/parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Celestial.dir/src/parser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/miquel/Desktop/Celestial/src/parser.cpp > CMakeFiles/Celestial.dir/src/parser.cpp.i

CMakeFiles/Celestial.dir/src/parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Celestial.dir/src/parser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/miquel/Desktop/Celestial/src/parser.cpp -o CMakeFiles/Celestial.dir/src/parser.cpp.s

# Object files for target Celestial
Celestial_OBJECTS = \
"CMakeFiles/Celestial.dir/src/main.cpp.o" \
"CMakeFiles/Celestial.dir/src/parser.cpp.o"

# External object files for target Celestial
Celestial_EXTERNAL_OBJECTS =

Celestial: CMakeFiles/Celestial.dir/src/main.cpp.o
Celestial: CMakeFiles/Celestial.dir/src/parser.cpp.o
Celestial: CMakeFiles/Celestial.dir/build.make
Celestial: CMakeFiles/Celestial.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/miquel/Desktop/Celestial/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable Celestial"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Celestial.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Celestial.dir/build: Celestial

.PHONY : CMakeFiles/Celestial.dir/build

CMakeFiles/Celestial.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Celestial.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Celestial.dir/clean

CMakeFiles/Celestial.dir/depend:
	cd /home/miquel/Desktop/Celestial/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/miquel/Desktop/Celestial /home/miquel/Desktop/Celestial /home/miquel/Desktop/Celestial/cmake-build-debug /home/miquel/Desktop/Celestial/cmake-build-debug /home/miquel/Desktop/Celestial/cmake-build-debug/CMakeFiles/Celestial.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Celestial.dir/depend

