cmake_minimum_required(VERSION 3.13.3)


# Create a header that includes all the ReadWrite States headers
# It also implements the ReadWriteMS constructor
#   This implementation adds all the previous formats/States to the machine.
#
# To be able to be automatically recognized, the new state should:
#   1.- Have the same Class Name as File Name (without the extension)
#   2.- Name must end with '_state'
#   3.- Be placed on ReadWriteMS folder

file(REMOVE "include/FormatHandling/formats_header.h")

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
        include/FormatHandling/*_state.h)

# Override formats_header,
#       add a header that informs this is auto-generated code
file(WRITE include/FormatHandling/formats_header.h
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
    file(APPEND include/FormatHandling/formats_header.h "#include \"${out}\"\n")

endforeach()

# Append the "#include" relative to the MachineStateHandler
# Append the begin of the constructor method
file(APPEND include/FormatHandling/formats_header.h
        "
#include \"FormatHandling/FormatManager.h\"

namespace FormatHandling {
FormatManager::FormatManager()
{
")

# Append to the file (to the constructor) a call
#       to add to the pool of states each state found
foreach(header_path ${stateheaders})
    get_filename_component(var1 ${header_path} NAME_WE)
    file(APPEND include/FormatHandling/formats_header.h "\taddState(new ${var1}(this));\n")
endforeach()

# Finish constructor method
file(APPEND include/FormatHandling/formats_header.h "};\n}")

# End header Creation
