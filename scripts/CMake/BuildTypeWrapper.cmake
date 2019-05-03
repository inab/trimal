cmake_minimum_required(VERSION 3.13.3)


set(default_build_type "Release")

# Extracted from https://blog.kitware.com/cmake-and-the-default-build-type/
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
else()
    message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}'.")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" lower_build_type)
if("${lower_build_type}" STREQUAL "debug")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wdouble-promotion -Wnull-dereference -Wimplicit-fallthrough=3 -Winline")
endif()