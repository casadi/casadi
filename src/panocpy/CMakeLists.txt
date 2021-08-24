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
            panoc-alm::casadi-loader
            pybind11::pybind11
    )
    target_compile_definitions(_panocpy PRIVATE VERSION_INFO=${PROJECT_VERSION})
    set(PANOCPY_RPATH "\${ORIGIN}/../../../.." "\${ORIGIN}")
    set_target_properties(_panocpy PROPERTIES
        DEBUG_POSTFIX ""
        ASAN_POSTFIX ""
        TSAN_POSTFIX ""
        CXX_VISIBILITY_PRESET hidden
        INSTALL_RPATH "${PANOCPY_RPATH}"
    )
    target_compile_definitions(_panocpy PRIVATE PANOCPY_MODULE_NAME=_panocpy)

    if (SKBUILD)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            execute_process(COMMAND bash -c "realpath \$(${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so)"
                            OUTPUT_VARIABLE libstdcpp_file
                            ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()
        install(TARGETS _panocpy DESTINATION .)
        install(CODE "set(libstdcpp_file \"${libstdcpp_file}\")")
        install(CODE [[
            file(GET_RUNTIME_DEPENDENCIES
                LIBRARIES $<TARGET_FILE:_panocpy>
                RESOLVED_DEPENDENCIES_VAR _r_deps
                UNRESOLVED_DEPENDENCIES_VAR _u_deps
                POST_EXCLUDE_REGEXES ".*"
                POST_INCLUDE_REGEXES "${CMAKE_BINARY_DIR}.*|.*libcasadi.*|.*libgcc.*"
            )
            foreach(_file ${_r_deps})
                file(INSTALL
                    DESTINATION "${CMAKE_INSTALL_PREFIX}"
                    TYPE SHARED_LIBRARY
                    FILES "${_file}"
                )
            endforeach()
            if(libstdcpp_file)
                file(INSTALL 
                    DESTINATION "${CMAKE_INSTALL_PREFIX}"
                    FILES "${libstdcpp_file}")
            endif()
            list(LENGTH _u_deps _u_length)
            if("${_u_length}" GREATER 0)
                message(WARNING "Unresolved dependencies detected!")
            endif()
        ]])
    endif()
endif()
