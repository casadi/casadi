find_package(Python3 REQUIRED COMPONENTS Interpreter)

macro(casadi_function_codegen_python target pythonfile)
    add_custom_command(OUTPUT "${target}.c"
                    COMMAND Python3::Interpreter 
                            "${CMAKE_CURRENT_SOURCE_DIR}/${pythonfile}"
                            "${target}"
                    MAIN_DEPENDENCY "${pythonfile}")
    add_library("${target}" SHARED "${target}.c")
    set_target_properties("${target}" PROPERTIES DEBUG_POSTFIX ""
                                                ASAN_POSTFIX ""
                                                TSAN_POSTFIX "")
    if (NOT WIN32)
        target_compile_options("${target}" PRIVATE "-Wno-unused-parameter")
    endif()
endmacro()
