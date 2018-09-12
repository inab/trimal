cmake_minimum_required(VERSION 2.8.12)
project(trimAl)
# Make project C++ 11
set(CMAKE_CXX_STANDARD 11)

# SET(BUILD_SHARED_LIBRARIES OFF)
# SET(CMAKE_EXE_LINKER_FLAGS "-static ")

# Make bin folder to output the compiled programs
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

#ADD_DEFINITIONS(-O2)

# Include the include folder. Yeah. I love redudancy
include_directories(include)

# Include all the cpp that are needed on our project
set(preview
        source/trimalManager.cpp
        source/newAlignment.cpp

        source/Statistics/statisticsMold.cpp
        source/Statistics/statisticsGaps.cpp
        source/Statistics/StatisticsManager.cpp
        source/Statistics/statisticsConservation.cpp
        source/Statistics/statisticsConsistency.cpp

        source/Cleaner.cpp
        source/ReadWriteMS/ReadWriteMachineState.cpp

        source/similarityMatrix.cpp
        source/sequencesMatrix.cpp

        source/utils.cpp

        source/reportsystem.cpp
        source/TimerFactory.cpp)

# Include all states of ReadWriteMachineState
FILE (GLOB states
        source/ReadWriteMS/*_state.cpp)

add_library(common OBJECT
        ${preview}
        ${states})

# Add the executables
add_executable(trimal
        source/newmain.cpp
        $<TARGET_OBJECTS:common>)

add_executable(readal
        source/ReadWriteMS/ReadWriteMSmain.cpp
        $<TARGET_OBJECTS:common>)

add_executable(sandbox
        source/sandbox.cpp
        $<TARGET_OBJECTS:common>)

# Link the mathematical library
target_link_libraries(trimal m)
target_link_libraries(readal m)

# Create a header that includes all the ReadWrite States headers
# It also implements the ReadWriteMS constructor
#   This implementation adds all the previous formats/States to the machine.
#
# To be able to be automatically recognized, the new state should:
#   1.- Have the same Class Name as File Name (without the extension)
#   2.- Name must end with '_state'
#   3.- Be placed on ReadWriteMS folder

message("Generating include/ReadWriteMS/formats_header.h related to format handling")

FILE (GLOB stateheaders
        include/ReadWriteMS/*_state.h)

file(WRITE include/ReadWriteMS/formats_header.h "//CMake generated code\n")

get_filename_component(var1 include ABSOLUTE)
foreach(header_path ${stateheaders})
    string(REPLACE "${var1}/" "" out "${header_path}")
    file(APPEND include/ReadWriteMS/formats_header.h "#include \"${out}\"\n")
endforeach()

file(APPEND include/ReadWriteMS/formats_header.h
"
#include \"ReadWriteMS/ReadWriteMachineState.h\"

ReadWriteMS::ReadWriteMS()
{
")

foreach(header_path ${stateheaders})
    get_filename_component(var1 ${header_path} NAME_WE)
    file(APPEND include/ReadWriteMS/formats_header.h "\taddState(new ${var1}(this));\n")
endforeach()

file(APPEND include/ReadWriteMS/formats_header.h "};")

# End header Creation
