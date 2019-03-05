cmake_minimum_required(VERSION 3.0.2)


# Specify the include folders, where all our headers are located.
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

# Create a lists of source files to be compiled into Object OBJLibraries
set(commonFiles
        source/Cleaner.cpp
        source/Alignment/Alignment.cpp
        source/Alignment/sequencesMatrix.cpp
        source/Statistics/similarityMatrix.cpp)

set(statisticFiles
        source/Statistics/Mold.cpp
        source/Statistics/Gaps.cpp
        source/Statistics/Manager.cpp
        source/Statistics/Similarity.cpp
        source/Statistics/Consistency.cpp)

set(reportSystemFiles
        source/reportsystem.cpp
        source/reportMessages/infoMessages.cpp
        source/reportMessages/errorMessages.cpp
        source/reportMessages/warningMessages.cpp)

FILE (GLOB FormatHandlerFiles
        source/FormatHandling/*_state.cpp)

FILE (GLOB tests
        tests/*.cpp
        tests/*/*.cpp)

# Create object libraries - This allows the system to reuse the same code for all targets
# They are splitted so they can be recompiled modularly an reduce recompilation time on change

# Core files of all projects, including the most basic code:
#   Alignment and dependencies
add_library(CoreOBJLib              OBJECT ${commonFiles})
# MSA Format handlers
add_library(FormatsOBJLib           OBJECT ${FormatHandlerFiles})
# MSA Format handler Manger
add_library(FormatHandlerOBJLib     OBJECT source/FormatHandling/BaseFormatHandler.cpp)
# Statistics files
add_library(StatisticOBJLib         OBJECT ${statisticFiles})
# Report system and error maps
add_library(ReportSystemOBJLib      OBJECT ${reportSystemFiles})
# Utils files
add_library(UtilsOBJLib             OBJECT source/utils.cpp)
# Internal Benchmarker
add_library(InternalBenchmarkOBJLib OBJECT source/InternalBenchmarker.cpp)

# Catch
add_library(CatchOBJLib             OBJECT tests/catch.hpp)
# Tests
add_library(TestsOBJLib             OBJECT ${tests})

