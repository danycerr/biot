# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /u/sw/pkgs/toolchains/gcc-glibc/5/base/bin/cmake

# The command to remove a file.
RM = /u/sw/pkgs/toolchains/gcc-glibc/5/base/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /u/archive/agip/cerroni/software/mygetfem/biot

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /u/archive/agip/cerroni/software/mygetfem/biot/build

# Include any dependencies generated for this target.
include CMakeFiles/biot.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/biot.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/biot.dir/flags.make

CMakeFiles/biot.dir/biot_ls.cpp.o: CMakeFiles/biot.dir/flags.make
CMakeFiles/biot.dir/biot_ls.cpp.o: ../biot_ls.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/archive/agip/cerroni/software/mygetfem/biot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/biot.dir/biot_ls.cpp.o"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/biot.dir/biot_ls.cpp.o -c /u/archive/agip/cerroni/software/mygetfem/biot/biot_ls.cpp

CMakeFiles/biot.dir/biot_ls.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/biot.dir/biot_ls.cpp.i"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/archive/agip/cerroni/software/mygetfem/biot/biot_ls.cpp > CMakeFiles/biot.dir/biot_ls.cpp.i

CMakeFiles/biot.dir/biot_ls.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/biot.dir/biot_ls.cpp.s"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/archive/agip/cerroni/software/mygetfem/biot/biot_ls.cpp -o CMakeFiles/biot.dir/biot_ls.cpp.s

CMakeFiles/biot.dir/biot_ls.cpp.o.requires:

.PHONY : CMakeFiles/biot.dir/biot_ls.cpp.o.requires

CMakeFiles/biot.dir/biot_ls.cpp.o.provides: CMakeFiles/biot.dir/biot_ls.cpp.o.requires
	$(MAKE) -f CMakeFiles/biot.dir/build.make CMakeFiles/biot.dir/biot_ls.cpp.o.provides.build
.PHONY : CMakeFiles/biot.dir/biot_ls.cpp.o.provides

CMakeFiles/biot.dir/biot_ls.cpp.o.provides.build: CMakeFiles/biot.dir/biot_ls.cpp.o


CMakeFiles/biot.dir/biot.cpp.o: CMakeFiles/biot.dir/flags.make
CMakeFiles/biot.dir/biot.cpp.o: ../biot.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/archive/agip/cerroni/software/mygetfem/biot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/biot.dir/biot.cpp.o"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/biot.dir/biot.cpp.o -c /u/archive/agip/cerroni/software/mygetfem/biot/biot.cpp

CMakeFiles/biot.dir/biot.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/biot.dir/biot.cpp.i"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/archive/agip/cerroni/software/mygetfem/biot/biot.cpp > CMakeFiles/biot.dir/biot.cpp.i

CMakeFiles/biot.dir/biot.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/biot.dir/biot.cpp.s"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/archive/agip/cerroni/software/mygetfem/biot/biot.cpp -o CMakeFiles/biot.dir/biot.cpp.s

CMakeFiles/biot.dir/biot.cpp.o.requires:

.PHONY : CMakeFiles/biot.dir/biot.cpp.o.requires

CMakeFiles/biot.dir/biot.cpp.o.provides: CMakeFiles/biot.dir/biot.cpp.o.requires
	$(MAKE) -f CMakeFiles/biot.dir/build.make CMakeFiles/biot.dir/biot.cpp.o.provides.build
.PHONY : CMakeFiles/biot.dir/biot.cpp.o.provides

CMakeFiles/biot.dir/biot.cpp.o.provides.build: CMakeFiles/biot.dir/biot.cpp.o


CMakeFiles/biot.dir/main.cpp.o: CMakeFiles/biot.dir/flags.make
CMakeFiles/biot.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/archive/agip/cerroni/software/mygetfem/biot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/biot.dir/main.cpp.o"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/biot.dir/main.cpp.o -c /u/archive/agip/cerroni/software/mygetfem/biot/main.cpp

CMakeFiles/biot.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/biot.dir/main.cpp.i"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/archive/agip/cerroni/software/mygetfem/biot/main.cpp > CMakeFiles/biot.dir/main.cpp.i

CMakeFiles/biot.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/biot.dir/main.cpp.s"
	/u/sw/pkgs/toolchains/gcc-glibc/5/prefix/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/archive/agip/cerroni/software/mygetfem/biot/main.cpp -o CMakeFiles/biot.dir/main.cpp.s

CMakeFiles/biot.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/biot.dir/main.cpp.o.requires

CMakeFiles/biot.dir/main.cpp.o.provides: CMakeFiles/biot.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/biot.dir/build.make CMakeFiles/biot.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/biot.dir/main.cpp.o.provides

CMakeFiles/biot.dir/main.cpp.o.provides.build: CMakeFiles/biot.dir/main.cpp.o


# Object files for target biot
biot_OBJECTS = \
"CMakeFiles/biot.dir/biot_ls.cpp.o" \
"CMakeFiles/biot.dir/biot.cpp.o" \
"CMakeFiles/biot.dir/main.cpp.o"

# External object files for target biot
biot_EXTERNAL_OBJECTS =

biot: CMakeFiles/biot.dir/biot_ls.cpp.o
biot: CMakeFiles/biot.dir/biot.cpp.o
biot: CMakeFiles/biot.dir/main.cpp.o
biot: CMakeFiles/biot.dir/build.make
biot: CMakeFiles/biot.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/u/archive/agip/cerroni/software/mygetfem/biot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable biot"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/biot.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/biot.dir/build: biot

.PHONY : CMakeFiles/biot.dir/build

CMakeFiles/biot.dir/requires: CMakeFiles/biot.dir/biot_ls.cpp.o.requires
CMakeFiles/biot.dir/requires: CMakeFiles/biot.dir/biot.cpp.o.requires
CMakeFiles/biot.dir/requires: CMakeFiles/biot.dir/main.cpp.o.requires

.PHONY : CMakeFiles/biot.dir/requires

CMakeFiles/biot.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/biot.dir/cmake_clean.cmake
.PHONY : CMakeFiles/biot.dir/clean

CMakeFiles/biot.dir/depend:
	cd /u/archive/agip/cerroni/software/mygetfem/biot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/archive/agip/cerroni/software/mygetfem/biot /u/archive/agip/cerroni/software/mygetfem/biot /u/archive/agip/cerroni/software/mygetfem/biot/build /u/archive/agip/cerroni/software/mygetfem/biot/build /u/archive/agip/cerroni/software/mygetfem/biot/build/CMakeFiles/biot.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/biot.dir/depend
