# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/daoce.wang/anaconda3/envs/cmake/bin/cmake

# The command to remove a file.
RM = /home/daoce.wang/anaconda3/envs/cmake/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3

# Include any dependencies generated for this target.
include test/CMakeFiles/convertBinToHDF5.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/convertBinToHDF5.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/convertBinToHDF5.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/convertBinToHDF5.dir/flags.make

test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o: test/CMakeFiles/convertBinToHDF5.dir/flags.make
test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o: test/src/convertBinToHDF5.cpp
test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o: test/CMakeFiles/convertBinToHDF5.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o"
	cd /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o -MF CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o.d -o CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o -c /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test/src/convertBinToHDF5.cpp

test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.i"
	cd /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test/src/convertBinToHDF5.cpp > CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.i

test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.s"
	cd /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test/src/convertBinToHDF5.cpp -o CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.s

# Object files for target convertBinToHDF5
convertBinToHDF5_OBJECTS = \
"CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o"

# External object files for target convertBinToHDF5
convertBinToHDF5_EXTERNAL_OBJECTS =

test/convertBinToHDF5: test/CMakeFiles/convertBinToHDF5.dir/src/convertBinToHDF5.o
test/convertBinToHDF5: test/CMakeFiles/convertBinToHDF5.dir/build.make
test/convertBinToHDF5: libhdf5sz3.a
test/convertBinToHDF5: /home/daoce.wang/hdf5-1.12.2/install/lib/libhdf5.so
test/convertBinToHDF5: test/CMakeFiles/convertBinToHDF5.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable convertBinToHDF5"
	cd /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/convertBinToHDF5.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/convertBinToHDF5.dir/build: test/convertBinToHDF5
.PHONY : test/CMakeFiles/convertBinToHDF5.dir/build

test/CMakeFiles/convertBinToHDF5.dir/clean:
	cd /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test && $(CMAKE_COMMAND) -P CMakeFiles/convertBinToHDF5.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/convertBinToHDF5.dir/clean

test/CMakeFiles/convertBinToHDF5.dir/depend:
	cd /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3 /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3 /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test /home/daoce.wang/develop/SZ_SLE/tools/H5Z-SZ3/test/CMakeFiles/convertBinToHDF5.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/convertBinToHDF5.dir/depend

