cmake_minimum_required(VERSION 2.8.12)
project(trimal)
# JEMALLOC
# Make project C++ 11
set(CMAKE_CXX_STANDARD 11)

# Make project libs static for use in Docker.
# Remember to comment these lines if you are building for non-docker use.
# SET(BUILD_SHARED_LIBRARIES OFF)
# SET(CMAKE_EXE_LINKER_FLAGS "-static")

# Make bin folder to output the compiled programs
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

ADD_DEFINITIONS(-O3)

# Include the include folder. Yeah. I love redudancy
include_directories(include)

# Include all the cpp that are needed on our project
set(preview
        source/trimalArgumentParser.cpp
        source/newAlignment.cpp

        source/Statistics/statisticsGaps.cpp
        source/Statistics/StatisticsManager.cpp
        source/Statistics/statisticsConservation2.cpp

        source/Cleaner.cpp
        source/ReadWriteMS/ReadWriteMachineState.cpp

        source/similarityMatrix.cpp
        source/sequencesMatrix.cpp

        source/utils.cpp
        source/compareFiles.cpp

        source/reportsystem.cpp
        source/TimerFactory.cpp
        )

# Include all states of ReadWriteMachineState
FILE (GLOB states
        source/ReadWriteMS/*_state.cpp)

FILE (GLOB test_files
        source/Catch/test_*.cpp)

SET (alignment_test_file source/Catch/test_generated_alignment.cpp)

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

add_executable(testMaker EXCLUDE_FROM_ALL
        source/Catch/TestMaker.cpp
        source/utils.cpp
        source/reportsystem.cpp)

add_executable(tests EXCLUDE_FROM_ALL
        source/Catch/main.cpp
        include/Catch/catch.hpp
        include/Catch/catchhelperfunctions.h
        source/Catch/Matchers/ArrayMatcher.cpp
        $<TARGET_OBJECTS:common>
        ${test_files}
        ${alignment_test_file})

# Command to build the alignment_test_files
add_dependencies(tests testMaker)

add_custom_command(OUTPUT ${alignment_test_file} COMMAND testMaker)

SET_TARGET_PROPERTIES(tests PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD True)

# Link the mathematical library
target_link_libraries(trimal m)
target_link_libraries(readal m)

# --------------------------- TESTING --------------------------- #

# find_program( MEMORYCHECK_COMMAND valgrind )
# set( MEMORYCHECK_COMMAND_OPTIONS "-v --tool=callgrind" )
# set( MEMORYCHECK_SUPPRESSIONS_FILE "${PROJECT_SOURCE_DIR}/valgrind_suppress.txt" )


# Basic tests
include(CTest)
add_test(Crash_NewTrimal    bin/trimal)
add_test(Crash_Readal       bin/readal)

# Macro that uses the compareOutput cmake to create tests
macro(compare_outputs test_name program_a program_b arguments)
    add_test(NAME ${test_name}
            COMMAND ${CMAKE_COMMAND}
            -Dprogram_a=${program_a}
            -Dprogram_b=${program_b}
            -Dcommand=${arguments}
            -P scripts/compareTrimAlOutput.cmake)
endmacro()

macro(compare_stdout test_name program_a program_b arguments)
    add_test(NAME ${test_name}
            COMMAND ${CMAKE_COMMAND}
            -Dprogram_a=${program_a}
            -Dprogram_b=${program_b}
            -Dcommand=${arguments}
            -P scripts/compareTrimAlSTD.cmake)
endmacro()

macro(compare_htmlsvg test_name program_a program_b arguments)
    add_test(NAME ${test_name}
            COMMAND ${CMAKE_COMMAND}
            -Dprogram_a=${program_a}
            -Dprogram_b=${program_b}
            -Dcommand=${arguments}
            -P scripts/compareTrimAlReports.cmake)
endmacro()

# Create lists of variables and files for testing
FILE(GLOB files
        dataset/*clw
        dataset/*phy
        dataset/*msl
        dataset/*fa*
        dataset/*pir
        dataset/*phy*
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
set(formats -nbrf -mega -nexus -clustal -fasta -fasta_m10 -phylip -phylip_m10 -phylip_paml -phylip_paml_m10 -phylip3.2 -phylip3.2_m10)
set(external -sgc -sgt -ssc -sst -sfc -sft -sident -soverlap)

set(VAR 3)

# Add all the combinations of arguments and files provided
foreach(file ${files})
    
    string(REGEX MATCH "[^/]+$" stripped_file ${file})
#     message("${stripped_file}")
    
    compare_outputs(
            "[.${VAR}.][trimal][file:${stripped_file}][vanilla]"
            "../TtrimAl/source/trimal" "bin/trimal"
            "-in ${file}")
    MATH(EXPR VAR "${VAR}+1")
            
    foreach(threshold ${thresholds})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][gt:${threshold}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} -gt ${threshold}")
        MATH(EXPR VAR "${VAR}+1")
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][st:${threshold}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} -st ${threshold}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()

    foreach(cons ${cons})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][cons:${cons}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} -cons ${cons}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()

    foreach(method ${automethod})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][${method}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} ${method}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()
    
    foreach(cluster ${clusters})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][cluster:${cluster}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} -clusters ${cluster}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()
    
    foreach(max_identity ${max_identities})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][max_identity:${max_identity}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} -maxidentity ${max_identity}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()
    
    foreach(block ${blocks})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][block:${block}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} -block ${block}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()

    foreach(seqoverlap ${seqoverlap})
        foreach(resoverlap ${resoverlap})
            compare_outputs(
                    "[.${VAR}.][trimal][file:${stripped_file}][seqo:${seqoverlap}][reso:${resoverlap}]"
                    "../TtrimAl/source/trimal" "bin/trimal"
                    "-in ${file} -resoverlap ${resoverlap} -seqoverlap ${seqoverlap}")
            MATH(EXPR VAR "${VAR}+1")
        endforeach()
    endforeach()
    
    foreach(format ${formats})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][format:${format}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} ${format}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()
    
    foreach(ext ${external})
        compare_outputs(
                "[.${VAR}.][trimal][file:${stripped_file}][externalParam:${ext}]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} ${ext}")
        MATH(EXPR VAR "${VAR}+1")
    endforeach()
    compare_htmlsvg(
                "[.${VAR}.][trimal][file:${stripped_file}][html]"
                "../TtrimAl/source/trimal" "bin/trimal"
                "-in ${file} -strict -htmlout")
        MATH(EXPR VAR "${VAR}+1")
        
endforeach()

compare_outputs(
            "[.${VAR}.][trimal][Compareset]"
            "../TtrimAl/source/trimal" "bin/trimal"
            "-compareset /home/vfernandez/git/trimal/dataset/alignments_comparison.X")
    MATH(EXPR VAR "${VAR}+1")
    
foreach(method ${automethod})
compare_outputs(
            "[.${VAR}.][trimal][Backtrans][${method}]"
            "../TtrimAl/source/trimal" "bin/trimal"
            "-in dataset/example.091.AA.strNOG.ENOG411BWBU.fasta -backtrans dataset/example.091.AA.strNOG.ENOG411BWBU.disordered.codon.fa ${method}")
    MATH(EXPR VAR "${VAR}+1")
endforeach()
    
    



# Force tests to be made each time new_trimal is built
# add_custom_command(
#        TARGET new_trimal
#       POST_BUILD
#      COMMAND ctest -R "vanilla")
