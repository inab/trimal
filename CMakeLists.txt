cmake_minimum_required(VERSION 2.8.12)
project(trimal)

# Make project C++ 11
set(CMAKE_CXX_STANDARD 11)

# Make bin folder to output the compiled programs
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# Set the debug flag
SET(CMAKE_BUILD_TYPE Debug)

# Include the include folder. Yeah. I love redudancy
include_directories(include)

# Include all the cpp that are needed on our project
set(preview
        source/alignment.cpp
        source/statisticsGaps.cpp
        source/utils.cpp
        source/similarityMatrix.cpp
        source/statisticsConservation.cpp
        source/sequencesMatrix.cpp
        source/compareFiles.cpp
        source/newAlignment.cpp
        source/Cleaner.cpp
        source/StatisticsManager.cpp
        source/ReadWriteManager.cpp
        source/trimalArgumentParser.cpp
        source/ReadWriteMS/ReadWriteMachineState.cpp
        )

# Include all states of ReadWriteMachineState
FILE (GLOB states source/ReadWriteMS/*_state.cpp)
# message("${states}")
# message("${preview}")

# add_library(states OBJECT ${states} ${preview})
add_library(common OBJECT ${preview} ${states})
        
# Add the executables
add_executable(trimal source/main.cpp $<TARGET_OBJECTS:common>)
add_executable(new_trimal source/newmain.cpp $<TARGET_OBJECTS:common>)
add_executable(readal source/readAl.cpp $<TARGET_OBJECTS:common>)
add_executable(statal source/statAl.cpp $<TARGET_OBJECTS:common>)
add_executable(readalMS source/ReadWriteMS/ReadWriteMSmain.cpp $<TARGET_OBJECTS:common>)

# Link the mathematical library
target_link_libraries(trimal m)
target_link_libraries(readal m)
target_link_libraries(statal m)
target_link_libraries(new_trimal m)
target_link_libraries(readalMS m)

# --------------------------- TESTING --------------------------- #

# new_trimal tests need trimal to be built
add_dependencies(new_trimal trimal)

# Basic tests
include(CTest)
add_test(Crash_Trimal       bin/trimal)
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
            -P Scripts/compareOutput.cmake)
endmacro()

# Create lists of variables and files for testing
FILE(GLOB files
        dataset/*clw
        dataset/*phy
       dataset/*msl
       dataset/*fa
       dataset/*pir
       dataset/*phy2
        )
set(thresholds 0 0.5 1)
set(cons 10 20 40 80)
set(automethod -gappyout -strict -strictplus -automated1 -nogaps -noallgaps)
set(resoverlap 0.1 0.3 0.6)
set(seqoverlap 5 10 50 80)
set(clusters 2 4)
set(max_identities 0 0.5 1)
set(blocks 2 4 6)

# Add all the combinations of gt and files provided
foreach(file ${files})
    
    string(REGEX MATCH "[^/]+$" stripped_file ${file})
#     message("${stripped_file}")
    
    compare_outputs(
            "[file:${stripped_file}]"
            "bin/trimal" "bin/new_trimal"
            "-in ${file}")
            
    foreach(threshold ${thresholds})
        compare_outputs(
                "[file:${stripped_file}][gt:${threshold}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -gt ${threshold}")
        compare_outputs(
                "[file:${stripped_file}][st:${threshold}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -st ${threshold}")
    endforeach()

    foreach(cons ${cons})
        compare_outputs(
                "[file:${stripped_file}][cons:${cons}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -cons ${cons}")
    endforeach()

    foreach(method ${automethod})
        compare_outputs(
                "[file:${stripped_file}][${method}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} ${method}")
    endforeach()
    
    foreach(cluster ${clusters})
        compare_outputs(
                "[file:${stripped_file}][cluster:${cluster}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -clusters ${cluster}")
    endforeach()
    
    foreach(max_identity ${max_identities})
        compare_outputs(
                "[file:${stripped_file}][max_identity:${max_identity}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -maxidentity ${max_identity}")
    endforeach()
    
    foreach(block ${blocks})
        compare_outputs(
                "[file:${stripped_file}][block:${block}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -block ${block}")
    endforeach()

    foreach(seqoverlap ${seqoverlap})
        foreach(resoverlap ${resoverlap})
            compare_outputs(
                    "[file:${stripped_file}][seqo:${seqoverlap}][reso:${resoverlap}]"
                    "bin/trimal" "bin/new_trimal"
                    "-in ${file} -resoverlap ${resoverlap} -seqoverlap ${seqoverlap}")
        endforeach()
    endforeach()
endforeach()

# # Force tests to be made each time new_trimal is built
# add_custom_command(
#        TARGET new_trimal
#       POST_BUILD
#      COMMAND ctest)
