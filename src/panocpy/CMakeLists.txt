if (SKBUILD)
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    execute_process(COMMAND "${Python3_EXECUTABLE}"
            -c "import pybind11; print(pybind11.get_include())"
        OUTPUT_VARIABLE panocpy_pybind11_include_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    add_library(pybind11::pybind11 INTERFACE IMPORTED)
    target_include_directories(pybind11::pybind11
        INTERFACE ${panocpy_pybind11_include_dir})

    Python3_add_library(_panocpy MODULE panocpy.cpp)
    target_link_libraries(_panocpy
        PRIVATE
            panoc-alm::panoc-alm
            pybind11::pybind11
    )
    target_compile_definitions(_panocpy PRIVATE VERSION_INFO=${PROJECT_VERSION})
    set_target_properties(_panocpy PROPERTIES
        DEBUG_POSTFIX ""
        ASAN_POSTFIX ""
        TSAN_POSTFIX ""
        CXX_VISIBILITY_PRESET hidden
    )
    target_compile_definitions(_panocpy PRIVATE PANOCPY_MODULE_NAME=_panocpy)

    install(TARGETS _panocpy DESTINATION .)
endif()
