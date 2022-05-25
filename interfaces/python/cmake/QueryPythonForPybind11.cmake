# First tries to find Python 3, then tries to import the pybind11 module to
# query the CMake config location, and finally imports pybind11 using
# find_package(pybind11 REQUIRED CONFIG).
function(find_pybind11_python_first)

    # Query Python to see if it knows where the headers are
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    if (NOT pybind11_ROOT OR NOT EXISTS ${pybind11_ROOT})
        execute_process(COMMAND ${Python3_EXECUTABLE}
                -m pybind11 --cmakedir
            OUTPUT_VARIABLE PY_BUILD_PYBIND11_CMAKE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE PY_BUILD_CMAKE_PYBIND11_RESULT)
        # If it was successful
        if (PY_BUILD_CMAKE_PYBIND11_RESULT EQUAL 0)
            message(STATUS "Found pybind11: ${PY_BUILD_PYBIND11_CMAKE}")
            set(pybind11_ROOT ${PY_BUILD_PYBIND11_CMAKE}
                CACHE PATH "Path to the pybind11 CMake configuration." FORCE)
        else()
            unset(pybind11_ROOT CACHE)
        endif()
    endif()

    find_package(pybind11 REQUIRED CONFIG)

endfunction()