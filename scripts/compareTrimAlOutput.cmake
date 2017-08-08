macro(compare program_a program_b command_list)
    
    file(REMOVE /tmp/in /tmp/out)
    
    execute_process(COMMAND (${program_a} ${command_list} -out /tmp/in) >> /tmp/null)
    execute_process(COMMAND (${program_b} ${command_list} -out /tmp/out) >> /tmp/null)
    
    if (NOT EXISTS /tmp/in AND EXISTS /tmp/out)
        message( SEND_ERROR "'${program_a}' didn't output a file, while '${program_b}' did." )
    elseif(EXISTS /tmp/in AND NOT EXISTS /tmp/out)
        message( SEND_ERROR "'${program_b}' didn't output a file, while '${program_a}' did." )
    elseif(NOT EXISTS /tmp/in AND NOT EXISTS /tmp/out)
        message( SEND_ERROR "'${program_b}' didn't output a file, neither '${program_a}'." )
    else()
        execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files /tmp/in /tmp/out RESULT_VARIABLE returnValue)
        if (${returnValue})
            message( SEND_ERROR "Result is not the same." )
        endif()
    endif()
    
endmacro(compare)
