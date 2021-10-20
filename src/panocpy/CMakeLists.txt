find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
execute_process(COMMAND "${Python3_EXECUTABLE}"
        -c "import pybind11; print(pybind11.get_include())"
    OUTPUT_VARIABLE panocpy_pybind11_include_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE PANOC_ALM_PYBIND11_RESULT)

if (PANOC_ALM_PYBIND11_RESULT EQUAL 0)
    add_library(pybind11::pybind11 INTERFACE IMPORTED)
    target_include_directories(pybind11::pybind11
        INTERFACE ${panocpy_pybind11_include_dir})

    Python3_add_library(_panocpy MODULE panocpy.cpp)
    target_link_libraries(_panocpy
        PRIVATE
            panoc-alm::panoc-alm-obj
            pybind11::pybind11
    )
    if (NOT WIN32)
        target_link_libraries(_panocpy 
            PRIVATE
                -static-libstdc++
                -static-libgcc
        )
    endif()
    target_compile_definitions(_panocpy PRIVATE VERSION_INFO=${PROJECT_VERSION})
    set(PANOCPY_RPATH "\${ORIGIN}/../../../.." "\${ORIGIN}")
    set_target_properties(_panocpy PROPERTIES
        DEBUG_POSTFIX ""
        ASAN_POSTFIX ""
        TSAN_POSTFIX ""
        CXX_VISIBILITY_PRESET hidden
        C_VISIBILITY_PRESET hidden
        INSTALL_RPATH "${PANOCPY_RPATH}"
    )
    if (NOT WIN32)
        target_link_options(_panocpy PRIVATE "LINKER:--exclude-libs,ALL")
    endif()
    target_compile_definitions(_panocpy PRIVATE PANOCPY_MODULE_NAME=_panocpy)
    if (TARGET panoc-alm::casadi-loader)
        target_compile_definitions(_panocpy PRIVATE PANOCPY_HAVE_CASADI)
        target_link_libraries(_panocpy PRIVATE panoc-alm::casadi-loader-obj)
    endif()

    option(PANOCPY_GEN_STUB "Generate Python stub files (.pyi) for the panocpy module." On)

    if (PANOCPY_GEN_STUB)
        set(STUBS_DIR ${CMAKE_CURRENT_BINARY_DIR}/stubs)
        add_custom_command(TARGET _panocpy POST_BUILD 
            COMMAND stubgen 
                    -m panocpy -m panocpy.casadi_problem -m panocpy._panocpy
                    -o ${STUBS_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
            USES_TERMINAL)
    endif()

    if (SKBUILD)
        install(TARGETS _panocpy
                DESTINATION .)
        if (PANOCPY_GEN_STUB)
            install(FILES ${STUBS_DIR}/panocpy/__init__.pyi
                          ${STUBS_DIR}/panocpy/casadi_problem.pyi
                          ${STUBS_DIR}/panocpy/_panocpy.pyi
                    DESTINATION .)
        endif()
    endif()
endif()
