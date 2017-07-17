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

# --------------------------- TESTING --------------------------- #

# new_trimal tests need trimal to be built
add_dependencies(new_trimal trimal)

# Basic tests
include(CTest)
add_test(Crash_Trimal bin/trimal)
add_test(Crash_Readal bin/readal)
add_test(Crash_Statal bin/statal)

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
set(gt 0 1)
set(st 0 1)
set(cons 60)
set(thold -gappyout -strict -automated1)
set(resoverlap 0.5)
set(seqoverlap 80)

# Add all the combinations of gt and files provided
foreach(file ${files})
    compare_outputs(
            "[file:${file}]"
            "bin/trimal" "bin/new_trimal"
            "-in ${file}")
    foreach(gt ${gt})
        compare_outputs(
                "[file:${file}][gt:${gt}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -gt ${gt}")
    endforeach()

    foreach(st ${st})
        compare_outputs(
                "[file:${file}][st:${st}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -st ${st}")
    endforeach()

    foreach(cons ${cons})
        compare_outputs(
                "[file:${file}][cons:${cons}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} -cons ${cons}")
    endforeach()

    foreach(thold ${thold})
        compare_outputs(
                "[file:${file}][${thold}]"
                "bin/trimal" "bin/new_trimal"
                "-in ${file} ${thold}")
    endforeach()

    foreach(seqoverlap ${seqoverlap})
        foreach(resoverlap ${resoverlap})
            compare_outputs(
                    "[file:${file}][seqo:${seqoverlap}][reso:${resoverlap}]"
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
