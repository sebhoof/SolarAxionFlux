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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.5/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.5/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build

# Include any dependencies generated for this target.
include src/CMakeFiles/test_axionflux.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/test_axionflux.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/test_axionflux.dir/flags.make

src/CMakeFiles/test_axionflux.dir/main.cpp.o: src/CMakeFiles/test_axionflux.dir/flags.make
src/CMakeFiles/test_axionflux.dir/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/test_axionflux.dir/main.cpp.o"
	cd /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/src && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_axionflux.dir/main.cpp.o -c /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/src/main.cpp

src/CMakeFiles/test_axionflux.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_axionflux.dir/main.cpp.i"
	cd /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/src/main.cpp > CMakeFiles/test_axionflux.dir/main.cpp.i

src/CMakeFiles/test_axionflux.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_axionflux.dir/main.cpp.s"
	cd /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/src && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/src/main.cpp -o CMakeFiles/test_axionflux.dir/main.cpp.s

# Object files for target test_axionflux
test_axionflux_OBJECTS = \
"CMakeFiles/test_axionflux.dir/main.cpp.o"

# External object files for target test_axionflux
test_axionflux_EXTERNAL_OBJECTS =

../bin/test_axionflux: src/CMakeFiles/test_axionflux.dir/main.cpp.o
../bin/test_axionflux: src/CMakeFiles/test_axionflux.dir/build.make
../bin/test_axionflux: ../lib/libaxionflux.so
../bin/test_axionflux: /usr/local/Cellar/gsl/2.6/lib/libgsl.dylib
../bin/test_axionflux: /usr/local/Cellar/gsl/2.6/lib/libgslcblas.dylib
../bin/test_axionflux: src/CMakeFiles/test_axionflux.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/test_axionflux"
	cd /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_axionflux.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/test_axionflux.dir/build: ../bin/test_axionflux

.PHONY : src/CMakeFiles/test_axionflux.dir/build

src/CMakeFiles/test_axionflux.dir/clean:
	cd /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/src && $(CMAKE_COMMAND) -P CMakeFiles/test_axionflux.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/test_axionflux.dir/clean

src/CMakeFiles/test_axionflux.dir/depend:
	cd /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/src /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/src /Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/build/src/CMakeFiles/test_axionflux.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/test_axionflux.dir/depend

