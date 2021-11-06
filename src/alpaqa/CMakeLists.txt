find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
execute_process(COMMAND "${Python3_EXECUTABLE}"
        -c "import pybind11; print(pybind11.get_include())"
    OUTPUT_VARIABLE alpaqa_pybind11_include_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE ALPAQA_PYBIND11_RESULT)

if (ALPAQA_PYBIND11_RESULT EQUAL 0)
    add_library(pybind11::pybind11 INTERFACE IMPORTED)
    target_include_directories(pybind11::pybind11
        INTERFACE ${alpaqa_pybind11_include_dir})

    Python3_add_library(_alpaqa MODULE alpaqa.cpp)
    target_link_libraries(_alpaqa
        PRIVATE
            alpaqa::alpaqa-obj
            pybind11::pybind11
    )
    if (NOT WIN32)
        target_link_libraries(_alpaqa 
            PRIVATE
                -static-libstdc++
                -static-libgcc
        )
    endif()
    target_compile_definitions(_alpaqa PRIVATE VERSION_INFO=${PROJECT_VERSION})
    set(ALPAQA_RPATH "\${ORIGIN}/../../../.." "\${ORIGIN}")
    set_target_properties(_alpaqa PROPERTIES
        DEBUG_POSTFIX ""
        ASAN_POSTFIX ""
        TSAN_POSTFIX ""
        CXX_VISIBILITY_PRESET hidden
        C_VISIBILITY_PRESET hidden
        INSTALL_RPATH "${ALPAQA_RPATH}"
    )
    if (NOT WIN32)
        target_link_options(_alpaqa PRIVATE "LINKER:--exclude-libs,ALL")
    endif()
    target_compile_definitions(_alpaqa PRIVATE ALPAQA_MODULE_NAME=_alpaqa)
    if (TARGET alpaqa::casadi-loader)
        target_compile_definitions(_alpaqa PRIVATE ALPAQA_HAVE_CASADI)
        target_link_libraries(_alpaqa PRIVATE alpaqa::casadi-loader-obj)
    endif()

    option(ALPAQA_GEN_STUB "Generate Python stub files (.pyi) for the alpaqa module." On)

    if (ALPAQA_GEN_STUB)
        find_program(STUBGEN_EXE stubgen REQUIRED)
        set(STUBS_DIR ${CMAKE_CURRENT_BINARY_DIR}/stubs)
        add_custom_command(TARGET _alpaqa POST_BUILD 
            COMMAND ${STUBGEN_EXE} 
                    -m alpaqa -m alpaqa.casadi_problem
                    -o ${STUBS_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
            USES_TERMINAL)
        add_custom_command(TARGET _alpaqa POST_BUILD 
            COMMAND ${STUBGEN_EXE} 
                    -m alpaqa._alpaqa
                    -o ${STUBS_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
            USES_TERMINAL)
    endif()

    if (SKBUILD)
        install(TARGETS _alpaqa
                DESTINATION .)
        if (ALPAQA_GEN_STUB)
            install(FILES ${STUBS_DIR}/alpaqa/__init__.pyi
                          ${STUBS_DIR}/alpaqa/casadi_problem.pyi
                          ${STUBS_DIR}/alpaqa/_alpaqa.pyi
                    DESTINATION .)
        endif()
    endif()
endif()
