macro(compare program_a program_b command_list)
    execute_process(COMMAND ${program_a} ${command_list} -out bin/in)
    execute_process(COMMAND ${program_b} ${command_list} -out bin/out)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files bin/in bin/out RESULT_VARIABLE returnValue)
    if (${returnValue})
        message( SEND_ERROR "Result is not the same." )
    endif()

#    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files bin/in.html bin/out.html RESULT_VARIABLE returnValue)
#    if (${returnValue})
#        message( SEND_ERROR "HTML result is not the same." )
#    endif()

    file(REMOVE bin/in.html bin/out.html bin/in bin/out)


endmacro(compare)

string(REPLACE " " ";" command_list ${command})
compare(${program_a} ${program_b} "${command_list}")