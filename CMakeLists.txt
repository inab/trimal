cmake_minimum_required(VERSION 3.13.3)
project(trimAl)

# Make project C++ 11
set(CMAKE_CXX_STANDARD 11)

# Specify bin folder to output the compiled programs
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# Set common files for all executables
set(TRIMAL_OBJECTS
    $<TARGET_OBJECTS:CoreOBJLib>
    $<TARGET_OBJECTS:FormatsOBJLib>
    $<TARGET_OBJECTS:FormatHandlerOBJLib>
    $<TARGET_OBJECTS:ReportSystemOBJLib>
    $<TARGET_OBJECTS:StatisticOBJLib>
    $<TARGET_OBJECTS:UtilsOBJLib>
    $<TARGET_OBJECTS:InternalBenchmarkOBJLib>)

# Scripts to detect various SIMD support for the current compiler
if (NOT DISABLE_AVX2)
    include("./scripts/CMake/FindAVX2.cmake")
endif()
if (NOT DISABLE_NEON)
    include("./scripts/CMake/FindNEON.cmake")
endif()
if (NOT DISABLE_SSE2)
    include("./scripts/CMake/FindSSE2.cmake")
endif()

# Add SIMD code if compiler supports it
if (HAVE_SSE2)
    list(APPEND TRIMAL_OBJECTS $<TARGET_OBJECTS:SSE2OBJLib>)
endif()
if (HAVE_AVX2)
    list(APPEND TRIMAL_OBJECTS $<TARGET_OBJECTS:AVX2OBJLib>)
endif()
if (HAVE_NEON)
    list(APPEND TRIMAL_OBJECTS $<TARGET_OBJECTS:NEONOBJLib>)
endif()

# Add the executables to the project.
add_executable(trimal
        source/trimAlMain.cpp               # Main
        source/trimalManager.cpp            # trimAl Manager

        source/VCFHandler.cpp               # VCF only present on trimAl
        ${TRIMAL_OBJECTS})

add_executable(readal
        source/FormatHandling/readAlMain.cpp        # Main

       ${TRIMAL_OBJECTS})

# Add the executables to the project.
add_executable(test
        tests/testMain.cpp
        tests/catch.hpp
        source/trimalManager.cpp            # trimAl Manager

        source/VCFHandler.cpp               # VCF only present on trimAl
        ${TRIMAL_OBJECTS})
SET_TARGET_PROPERTIES(test PROPERTIES EXCLUDE_FROM_ALL True)

# Add `cpu_features` for detecting supported SIMD at compile time
add_subdirectory(vendor/cpu_features)
include_directories(vendor/cpu_features/include)
target_link_libraries(trimal cpu_features)
target_link_libraries(readal cpu_features)

# Link the mathematical library to the targets
target_link_libraries(trimal m)
target_link_libraries(readal m)

target_link_libraries(test stdc++fs)

# Script that sets Release type to default and prints information on configuration
include("./scripts/CMake/BuildTypeWrapper.cmake")

# Script to handle static and dynamic compilation.
#   Static compilation is useful if compiling for containers without SO
include("./scripts/CMake/staticCompilationWrapper.cmake")

# Script to create intermediate OBJLib objects to include in targets
include("./scripts/CMake/OBJ-LIB-creator.cmake")

# Script to configure the overwrite policy of the format manager:
#   - Overwrite and warn.
#   - Dont overwrite unless needed. If needed, overwrite original.
#   - Dont overwrite unless needed. If needed, overwrite last available.
include("./scripts/CMake/FormatHandlerOverwritePolicy.cmake")

# Script to create the Format Handler constructor, which incorporates
#   all found stats to the internal pool of the handler
include("./scripts/CMake/FormatHandlerHeaderCreator.cmake")

