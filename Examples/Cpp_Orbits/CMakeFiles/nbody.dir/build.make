# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits

# Include any dependencies generated for this target.
include CMakeFiles/nbody.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/nbody.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/nbody.dir/flags.make

CMakeFiles/nbody.dir/simulation.cpp.o: CMakeFiles/nbody.dir/flags.make
CMakeFiles/nbody.dir/simulation.cpp.o: simulation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/nbody.dir/simulation.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nbody.dir/simulation.cpp.o -c /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/simulation.cpp

CMakeFiles/nbody.dir/simulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nbody.dir/simulation.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/simulation.cpp > CMakeFiles/nbody.dir/simulation.cpp.i

CMakeFiles/nbody.dir/simulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nbody.dir/simulation.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/simulation.cpp -o CMakeFiles/nbody.dir/simulation.cpp.s

CMakeFiles/nbody.dir/simulation.cpp.o.requires:

.PHONY : CMakeFiles/nbody.dir/simulation.cpp.o.requires

CMakeFiles/nbody.dir/simulation.cpp.o.provides: CMakeFiles/nbody.dir/simulation.cpp.o.requires
	$(MAKE) -f CMakeFiles/nbody.dir/build.make CMakeFiles/nbody.dir/simulation.cpp.o.provides.build
.PHONY : CMakeFiles/nbody.dir/simulation.cpp.o.provides

CMakeFiles/nbody.dir/simulation.cpp.o.provides.build: CMakeFiles/nbody.dir/simulation.cpp.o


CMakeFiles/nbody.dir/orbit_integration.cpp.o: CMakeFiles/nbody.dir/flags.make
CMakeFiles/nbody.dir/orbit_integration.cpp.o: orbit_integration.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/nbody.dir/orbit_integration.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nbody.dir/orbit_integration.cpp.o -c /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/orbit_integration.cpp

CMakeFiles/nbody.dir/orbit_integration.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nbody.dir/orbit_integration.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/orbit_integration.cpp > CMakeFiles/nbody.dir/orbit_integration.cpp.i

CMakeFiles/nbody.dir/orbit_integration.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nbody.dir/orbit_integration.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/orbit_integration.cpp -o CMakeFiles/nbody.dir/orbit_integration.cpp.s

CMakeFiles/nbody.dir/orbit_integration.cpp.o.requires:

.PHONY : CMakeFiles/nbody.dir/orbit_integration.cpp.o.requires

CMakeFiles/nbody.dir/orbit_integration.cpp.o.provides: CMakeFiles/nbody.dir/orbit_integration.cpp.o.requires
	$(MAKE) -f CMakeFiles/nbody.dir/build.make CMakeFiles/nbody.dir/orbit_integration.cpp.o.provides.build
.PHONY : CMakeFiles/nbody.dir/orbit_integration.cpp.o.provides

CMakeFiles/nbody.dir/orbit_integration.cpp.o.provides.build: CMakeFiles/nbody.dir/orbit_integration.cpp.o


# Object files for target nbody
nbody_OBJECTS = \
"CMakeFiles/nbody.dir/simulation.cpp.o" \
"CMakeFiles/nbody.dir/orbit_integration.cpp.o"

# External object files for target nbody
nbody_EXTERNAL_OBJECTS =

nbody: CMakeFiles/nbody.dir/simulation.cpp.o
nbody: CMakeFiles/nbody.dir/orbit_integration.cpp.o
nbody: CMakeFiles/nbody.dir/build.make
nbody: CMakeFiles/nbody.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable nbody"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nbody.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/nbody.dir/build: nbody

.PHONY : CMakeFiles/nbody.dir/build

CMakeFiles/nbody.dir/requires: CMakeFiles/nbody.dir/simulation.cpp.o.requires
CMakeFiles/nbody.dir/requires: CMakeFiles/nbody.dir/orbit_integration.cpp.o.requires

.PHONY : CMakeFiles/nbody.dir/requires

CMakeFiles/nbody.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/nbody.dir/cmake_clean.cmake
.PHONY : CMakeFiles/nbody.dir/clean

CMakeFiles/nbody.dir/depend:
	cd /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits /home/miquel/Desktop/Celestial/Examples/Cpp_Orbits/CMakeFiles/nbody.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/nbody.dir/depend
