
macro(compare program_a program_b command_list)

 string(RANDOM RNUMBER LENGTH 20)
#  message(${RNUMBER})

    execute_process(COMMAND time ${program_a} ${command_list} -out "/tmp/in${RNUMBER}" )
    
    if (NOT EXISTS /tmp/in)
        message ( WARNING "Original algorithm failed to make an output file. This may be caused by an empty alignment after cleanse. Don't trying with updated algorithm." )
    
    else()
        execute_process(COMMAND ${program_b} ${command_list} -out "/tmp/out${RNUMBER}" )
        execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "/tmp/in${RNUMBER}" "/tmp/out${RNUMBER}" RESULT_VARIABLE returnValue)
        
        if (${returnValue})
            message( SEND_ERROR "Result is not the same." )
        endif()
    endif()
    
#     file(REMOVE "tmp/in${RNUMBER}" "tmp/out${RNUMBER}")
endmacro(compare)

string(REPLACE " " ";" command_list ${command})
compare(${program_a} ${program_b} "${command_list}")
