
macro(compare program_a program_b command_list)

 string(RANDOM LENGTH 20 RNUMBER)

    execute_process(COMMAND ${program_a} ${command_list} -out "/tmp/in${RNUMBER}" )
    
    if (NOT EXISTS /tmp/in${RNUMBER})
        message ( WARNING "Original algorithm failed to make an output file. This may be caused by an empty alignment after cleanse. Don't trying with updated algorithm." )

    else()
        execute_process(COMMAND ${program_b} ${command_list} -out "/tmp/out${RNUMBER}" )
        execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "/tmp/in${RNUMBER}" "/tmp/out${RNUMBER}" RESULT_VARIABLE returnValue)
        
        if (${returnValue})
            message( SEND_ERROR "Result is not the same. ${RNUMBER}" )
            message( SEND_ERROR "meld /tmp/in${RNUMBER} /tmp/out${RNUMBER} &" )
        else()
            file(REMOVE "tmp/in${RNUMBER}" "tmp/out${RNUMBER}")
        endif()
        
    endif()
    
    
endmacro(compare)

string(REPLACE " " ";" command_list ${command})
compare(${program_a} ${program_b} "${command_list}")

