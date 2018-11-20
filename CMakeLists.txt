cmake_minimum_required(VERSION 2.8.12)
project(trimAl)

# Make project C++ 11
set(CMAKE_CXX_STANDARD 11)

# This is needed to create a static binary,
#       which can be used directly into a docker
#       without adding any other dependency.
# Not recommended if not creating a docker image.

#SET(BUILD_SHARED_LIBRARIES OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")




# Specify bin folder to output the compiled programs
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# Specify the level of optimization used
#ADD_DEFINITIONS(-O2)




# Specify the include folders, where all our headers are located.
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")


# Create a list of all common source files
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
        source/timerFactory.cpp)

# Create a list of all format states source files
FILE (GLOB states
        source/ReadWriteMS/*_state.cpp)

# Create a library including all states and common source files
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

message(
"-- Generating include/ReadWriteMS/formats_header.h related to format handling.
   This file includes all States found, and also defines the ReadWriteMS constructor.
   To be able to be automatically recognized, the new state should:
        1.- Have the same Class Name as File Name (without the extension)
        2.- Name must end with '_state'
        3.- Be placed on ReadWriteMS folder"
)

#  Create a list of all format states header files
FILE (GLOB stateheaders
        include/ReadWriteMS/*_state.h)

# Override formats_header, 
#       add a header that informs this is auto-generated code
file(WRITE include/ReadWriteMS/formats_header.h 
"// CMake generated code
// Do not manually modify this file

// This file includes all States found, and also defines the ReadWriteMS constructor.
// To be able to be automatically recognized, the new state should:
//        1.- Have the same Class Name as File Name (without the extension)
//        2.- Name must end with '_state'
//        3.- Be placed on ReadWriteMS folder

")

# Get the absolute path to include folder
get_filename_component(var1 include ABSOLUTE)


foreach(header_path ${stateheaders})

    # Remove absolute path to include folders 
    string(REPLACE "${var1}/" "" out "${header_path}")

    # Append the "#include X" where X is the route to the header, 
    #   not including the absolute path to include folder
    file(APPEND include/ReadWriteMS/formats_header.h "#include \"${out}\"\n")

endforeach()

# Append the "#include" relative to the MachineStateHandler
# Append the begin of the constructor method
file(APPEND include/ReadWriteMS/formats_header.h
"
#include \"ReadWriteMS/ReadWriteMachineState.h\"

ReadWriteMS::ReadWriteMS()
{
")

# Append to the file (to the constructor) a call 
#       to add to the pool of states each state found
foreach(header_path ${stateheaders})
    get_filename_component(var1 ${header_path} NAME_WE)
    file(APPEND include/ReadWriteMS/formats_header.h "\taddState(new ${var1}(this));\n")
endforeach()

# Finish constructor method
file(APPEND include/ReadWriteMS/formats_header.h "};")

# End header Creation
