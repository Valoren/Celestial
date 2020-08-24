# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /snap/clion/123/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/123/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/miquel/Desktop/Celestial

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/miquel/Desktop/Celestial

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/snap/clion/123/bin/cmake/linux/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/snap/clion/123/bin/cmake/linux/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/miquel/Desktop/Celestial/CMakeFiles /home/miquel/Desktop/Celestial/CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/miquel/Desktop/Celestial/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named Celestial

# Build rule for target.
Celestial: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 Celestial
.PHONY : Celestial

# fast build rule for target.
Celestial/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/build
.PHONY : Celestial/fast

src/integration.o: src/integration.cpp.o

.PHONY : src/integration.o

# target to build an object file
src/integration.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/integration.cpp.o
.PHONY : src/integration.cpp.o

src/integration.i: src/integration.cpp.i

.PHONY : src/integration.i

# target to preprocess a source file
src/integration.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/integration.cpp.i
.PHONY : src/integration.cpp.i

src/integration.s: src/integration.cpp.s

.PHONY : src/integration.s

# target to generate assembly for a file
src/integration.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/integration.cpp.s
.PHONY : src/integration.cpp.s

src/main.o: src/main.cpp.o

.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i

.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s

.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

src/menu.o: src/menu.cpp.o

.PHONY : src/menu.o

# target to build an object file
src/menu.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/menu.cpp.o
.PHONY : src/menu.cpp.o

src/menu.i: src/menu.cpp.i

.PHONY : src/menu.i

# target to preprocess a source file
src/menu.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/menu.cpp.i
.PHONY : src/menu.cpp.i

src/menu.s: src/menu.cpp.s

.PHONY : src/menu.s

# target to generate assembly for a file
src/menu.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/menu.cpp.s
.PHONY : src/menu.cpp.s

src/parser.o: src/parser.cpp.o

.PHONY : src/parser.o

# target to build an object file
src/parser.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/parser.cpp.o
.PHONY : src/parser.cpp.o

src/parser.i: src/parser.cpp.i

.PHONY : src/parser.i

# target to preprocess a source file
src/parser.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/parser.cpp.i
.PHONY : src/parser.cpp.i

src/parser.s: src/parser.cpp.s

.PHONY : src/parser.s

# target to generate assembly for a file
src/parser.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/parser.cpp.s
.PHONY : src/parser.cpp.s

src/vector.o: src/vector.cpp.o

.PHONY : src/vector.o

# target to build an object file
src/vector.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/vector.cpp.o
.PHONY : src/vector.cpp.o

src/vector.i: src/vector.cpp.i

.PHONY : src/vector.i

# target to preprocess a source file
src/vector.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/vector.cpp.i
.PHONY : src/vector.cpp.i

src/vector.s: src/vector.cpp.s

.PHONY : src/vector.s

# target to generate assembly for a file
src/vector.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Celestial.dir/build.make CMakeFiles/Celestial.dir/src/vector.cpp.s
.PHONY : src/vector.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... Celestial"
	@echo "... src/integration.o"
	@echo "... src/integration.i"
	@echo "... src/integration.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/menu.o"
	@echo "... src/menu.i"
	@echo "... src/menu.s"
	@echo "... src/parser.o"
	@echo "... src/parser.i"
	@echo "... src/parser.s"
	@echo "... src/vector.o"
	@echo "... src/vector.i"
	@echo "... src/vector.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

