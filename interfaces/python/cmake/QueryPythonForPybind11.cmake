option(USE_GLOBAL_PYBIND11 "Don't query Python to find pybind11" Off)
mark_as_advanced(USE_GLOBAL_PYBIND11)

# First tries to find Python 3, then tries to import the pybind11 module to
# query the CMake config location, and finally imports pybind11 using
# find_package(pybind11 REQUIRED CONFIG CMAKE_FIND_ROOT_PATH_BOTH).
function(find_pybind11_python_first)

    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    if (NOT USE_GLOBAL_PYBIND11)
        # Query Python to see if it knows where the headers are
        if (NOT pybind11_ROOT OR NOT EXISTS ${pybind11_ROOT})
            message(STATUS "Detecting pybind11 CMake location")
            execute_process(COMMAND ${Python3_EXECUTABLE}
                    -m pybind11 --cmakedir
                OUTPUT_VARIABLE PY_BUILD_PYBIND11_ROOT
                OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE PY_BUILD_CMAKE_PYBIND11_RESULT)
            # If it was successful
            if (PY_BUILD_CMAKE_PYBIND11_RESULT EQUAL 0)
                message(STATUS "pybind11 CMake location: ${PY_BUILD_PYBIND11_ROOT}")
                set(pybind11_ROOT ${PY_BUILD_PYBIND11_ROOT}
                    CACHE PATH "Path to the pybind11 CMake configuration." FORCE)
            else()
                unset(pybind11_ROOT CACHE)
            endif()
        endif()
    endif()

    # pybind11 is header-only, so finding a native version is fine
    find_package(pybind11 ${ARGN} REQUIRED CONFIG CMAKE_FIND_ROOT_PATH_BOTH)

    # Tweak extension suffix when cross-compiling
    if (CMAKE_CROSSCOMPILING)
        if (NOT PY_BUILD_EXT_SUFFIX)
            message(STATUS "Determining Python extension suffix")
            # Find the python3.x-config script in the sysroot instead of on the
            # build system:
            find_program(PY_BUILD_Python3_CONFIG
                python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}-config
                ONLY_CMAKE_FIND_ROOT_PATH)
            # Report errors:
            if (NOT PY_BUILD_Python3_CONFIG)
                message(FATAL_ERROR "Unable to find python3-config."
                    "\nTry manually setting PY_BUILD_EXT_SUFFIX.")
            else()
                # If we found the python3.x-config script, query it for the
                # extension suffix:
                execute_process(COMMAND ${PY_BUILD_Python3_CONFIG}
                    --extension-suffix
                    OUTPUT_VARIABLE PY_BUILD_EXT_SUFFIX
                    ERROR_VARIABLE PY_BUILD_EXT_SUFFIX_ERR
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE PY_BUILD_EXT_SUFFIX_RESULT)
                # Report errors:
                if (NOT PY_BUILD_EXT_SUFFIX_RESULT EQUAL 0
                    OR NOT PY_BUILD_EXT_SUFFIX)
                    message(FATAL_ERROR "Unable to determine extension suffix:"
                        "\n${PY_BUILD_EXT_SUFFIX}"
                        "\n${PY_BUILD_EXT_SUFFIX_ERR}"
                        "\nTry manually setting PY_BUILD_EXT_SUFFIX.")
                endif()
                # Cache the result:
                set(PY_BUILD_EXT_SUFFIX ${PY_BUILD_EXT_SUFFIX} CACHE STRING
                    "The extension for Python extension modules")
            endif()
        endif()
        # Override pybind11NewTools.cmake's PYTHON_MODULE_EXTENSION variable:
        message(STATUS "Python extension suffix: ${PY_BUILD_EXT_SUFFIX}")
        set(PYTHON_MODULE_EXTENSION ${PY_BUILD_EXT_SUFFIX}
            CACHE INTERNAL "" FORCE)
    endif()

endfunction()