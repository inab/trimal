cmake_minimum_required(VERSION 2.8.12)
project(trimal)

# Make project C++ 11
set(CMAKE_CXX_STANDARD 11)

# Make project libs static for use in Docker.
# Remember to comment these lines if you are building for non-docker use.
# SET(BUILD_SHARED_LIBRARIES OFF)
# SET(CMAKE_EXE_LINKER_FLAGS "-static")

# Make bin folder to output the compiled programs
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# Set the debug flag
SET(CMAKE_BUILD_TYPE Debug)

# ADD_DEFINITIONS(-g -O2 -fsigned-char -freg-struct-return -Wall -W -Wshadow -Wstrict-prototypes -Wpointer-arith -Wcast-qual -Winline )
ADD_DEFINITIONS(-Winline -march=native --param max-inline-insns-single=4000 --param inline-unit-growth=200 )

# Include the include folder. Yeah. I love redudancy
include_directories(include)

# Include all the cpp that are needed on our project
set(preview
        source/statisticsGaps.cpp
        source/utils.cpp
        source/similarityMatrix.cpp
        source/statisticsConservation.cpp
        source/sequencesMatrix.cpp
        source/compareFiles.cpp
        source/newAlignment.cpp
        source/Cleaner.cpp
        source/StatisticsManager.cpp
        source/trimalArgumentParser.cpp
        source/ReadWriteMS/ReadWriteMachineState.cpp
        source/verbosemanager.cpp
        )

# Include all states of ReadWriteMachineState
FILE (GLOB states source/ReadWriteMS/*_state.cpp)

# add_library(states OBJECT ${states} ${preview})
add_library(common OBJECT ${preview} ${states})
        
# Add the executables
# add_executable(trimal source/main.cpp $<TARGET_OBJECTS:common>)
add_executable(trimal source/newmain.cpp $<TARGET_OBJECTS:common>)
# add_executable(readal source/readAl.cpp $<TARGET_OBJECTS:common>)
# add_executable(statal source/statAl.cpp $<TARGET_OBJECTS:common>)
add_executable(readal source/ReadWriteMS/ReadWriteMSmain.cpp $<TARGET_OBJECTS:common>)

# Link the mathematical library
target_link_libraries(trimal m)
target_link_libraries(readal m)
# target_link_libraries(statal m)
# target_link_libraries(new_trimal m)
# target_link_libraries(readalMS m)

# --------------------------- TESTING --------------------------- #

find_program( MEMORYCHECK_COMMAND valgrind )
set( MEMORYCHECK_COMMAND_OPTIONS "-v --tool=callgrind" )
# set( MEMORYCHECK_SUPPRESSIONS_FILE "${PROJECT_SOURCE_DIR}/valgrind_suppress.txt" )


# Basic tests
include(CTest)
# set_tests_properties(PROPERTIES TIMEOUT 0)
add_test(Crash_Trimal       ../Strimal/source/trimal)
add_test(Crash_NewTrimal    bin/new_trimal)
add_test(Crash_Readal       bin/readal)
add_test(Crash_Statal       bin/statal)
add_test(Crash_Readal       bin/readalMS)

# Macro that uses the compareOutput cmake to create tests
macro(compare_outputs test_name program_a program_b arguments)
    add_test(NAME ${test_name}
            COMMAND ${CMAKE_COMMAND}
            -Dprogram_a=${program_a}
            -Dprogram_b=${program_b}
            -Dcommand=${arguments}
            -P scripts/compareTrimAlOutput.cmake)
endmacro()

# Create lists of variables and files for testing
FILE(GLOB files
        dataset/*clw
        dataset/*phy
        dataset/*msl
        dataset/*fa*
        dataset/*pir
        dataset/*phy2
    )
    
# Create a pack of arguments to use for testing
set(thresholds 0.5)
set(cons 10)
set(automethod -gappyout -strict -strictplus -automated1 -nogaps -noallgaps)
set(resoverlap 0.5)
set(seqoverlap 5)
set(clusters 2)
set(max_identities 0.5)
set(blocks 2)

# Add all the combinations of arguments and files provided
foreach(file ${files})
    
    string(REGEX MATCH "[^/]+$" stripped_file ${file})
#     message("${stripped_file}")
    
    compare_outputs(
            "[trimal][file:${stripped_file}][vanilla]"
            "../Strimal/source/trimal" "bin/trimal"
            "-in ${file}")
            
    foreach(threshold ${thresholds})
        compare_outputs(
                "[trimal][file:${stripped_file}][gt:${threshold}]"
                "../Strimal/source/trimal" "bin/trimal"
                "-in ${file} -gt ${threshold}")
        compare_outputs(
                "[trimal][file:${stripped_file}][st:${threshold}]"
                "../Strimal/source/trimal" "bin/trimal"
                "-in ${file} -st ${threshold}")
    endforeach()

    foreach(cons ${cons})
        compare_outputs(
                "[trimal][file:${stripped_file}][cons:${cons}]"
                "../Strimal/source/trimal" "bin/trimal"
                "-in ${file} -cons ${cons}")
    endforeach()

    foreach(method ${automethod})
        compare_outputs(
                "[trimal][file:${stripped_file}][${method}]"
                "../Strimal/source/trimal" "bin/trimal"
                "-in ${file} ${method}")
    endforeach()
    
    foreach(cluster ${clusters})
        compare_outputs(
                "[trimal][file:${stripped_file}][cluster:${cluster}]"
                "../Strimal/source/trimal" "bin/trimal"
                "-in ${file} -clusters ${cluster}")
    endforeach()
    
    foreach(max_identity ${max_identities})
        compare_outputs(
                "[trimal][file:${stripped_file}][max_identity:${max_identity}]"
                "../Strimal/source/trimal" "bin/trimal"
                "-in ${file} -maxidentity ${max_identity}")
    endforeach()
    
    foreach(block ${blocks})
        compare_outputs(
                "[trimal][file:${stripped_file}][block:${block}]"
                "../Strimal/source/trimal" "bin/trimal"
                "-in ${file} -block ${block}")
    endforeach()

    foreach(seqoverlap ${seqoverlap})
        foreach(resoverlap ${resoverlap})
            compare_outputs(
                    "[trimal][file:${stripped_file}][seqo:${seqoverlap}][reso:${resoverlap}]"
                    "../Strimal/source/trimal" "bin/trimal"
                    "-in ${file} -resoverlap ${resoverlap} -seqoverlap ${seqoverlap}")
        endforeach()
    endforeach()
endforeach()

# Force tests to be made each time new_trimal is built
# add_custom_command(
#        TARGET new_trimal
#       POST_BUILD
#      COMMAND ctest -R "vanilla")
