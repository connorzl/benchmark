# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.12.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.12.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Connor/research/examples/basic_project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Connor/research/examples/basic_project/build

# Include any dependencies generated for this target.
include deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/depend.make

# Include the progress variables for this target.
include deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/progress.make

# Include the compile flags for this target's objects.
include deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/flags.make

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui.cpp.o: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/flags.make
deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui.cpp.o: ../deps/polyscope/deps/imgui/imgui/imgui.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Connor/research/examples/basic_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui.cpp.o"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui/imgui.cpp.o -c /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui.cpp

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui/imgui.cpp.i"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui.cpp > CMakeFiles/imgui.dir/imgui/imgui.cpp.i

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui/imgui.cpp.s"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui.cpp -o CMakeFiles/imgui.dir/imgui/imgui.cpp.s

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.o: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/flags.make
deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.o: ../deps/polyscope/deps/imgui/imgui/imgui_draw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Connor/research/examples/basic_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.o"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.o -c /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui_draw.cpp

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.i"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui_draw.cpp > CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.i

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.s"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui_draw.cpp -o CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.s

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.o: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/flags.make
deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.o: ../deps/polyscope/deps/imgui/imgui/imgui_demo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Connor/research/examples/basic_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.o"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.o -c /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui_demo.cpp

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.i"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui_demo.cpp > CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.i

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.s"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui/imgui/imgui_demo.cpp -o CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.s

# Object files for target imgui
imgui_OBJECTS = \
"CMakeFiles/imgui.dir/imgui/imgui.cpp.o" \
"CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.o" \
"CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.o"

# External object files for target imgui
imgui_EXTERNAL_OBJECTS =

lib/libimgui.a: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui.cpp.o
lib/libimgui.a: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_draw.cpp.o
lib/libimgui.a: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/imgui/imgui_demo.cpp.o
lib/libimgui.a: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/build.make
lib/libimgui.a: deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Connor/research/examples/basic_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library ../../../../lib/libimgui.a"
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean_target.cmake
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/imgui.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/build: lib/libimgui.a

.PHONY : deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/build

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/clean:
	cd /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean.cmake
.PHONY : deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/clean

deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/depend:
	cd /Users/Connor/research/examples/basic_project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Connor/research/examples/basic_project /Users/Connor/research/examples/basic_project/deps/polyscope/deps/imgui /Users/Connor/research/examples/basic_project/build /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui /Users/Connor/research/examples/basic_project/build/deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : deps/polyscope/deps/imgui/CMakeFiles/imgui.dir/depend

