cmake_minimum_required(VERSION 3.13.3)


# This is needed to create a static binary,
#       which can be used directly into a docker
#       without adding any other dependency.
# Not recommended if not creating a docker image.
#SET(BUILD_SHARED_LIBRARIES OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

set(static false CACHE BOOL "Compile the targets statically? This allows to run the targets on containers without including a SO.")
if (${static})
    message(STATUS "Configuring targets to compile 'Statically'.")
    set(BUILD_SHARED_LIBRARIES OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
else()
    message(STATUS "Configuring targets to compile 'Dynamically'.")
endif()