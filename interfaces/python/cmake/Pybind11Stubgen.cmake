cmake_minimum_required(VERSION 4.1)

function(pybind11_stubgen target)

    # Parse arguments
    set(OPTIONS
        # Package folder containing the Python extension module, relative to the
        # installation prefix (CMAKE_INSTALL_PREFIX). This should match the
        # DESTINATION argument of the given target.
        PACKAGE
        # The CMake installation component that the stub generation should be
        # part of.
        COMPONENT
        # Relative path of the Python package in the installation prefix.
        PACKAGE_ROOT
    )
    cmake_parse_arguments(STUBGEN "" "${OPTIONS}" "" ${ARGN})
    if (NOT DEFINED STUBGEN_PACKAGE)
        set(STUBGEN_PACKAGE ${PY_BUILD_CMAKE_IMPORT_NAME})
    endif()
    if (NOT DEFINED STUBGEN_COMPONENT)
        set(STUBGEN_COMPONENT "python_stubs")
    endif()
    if (NOT DEFINED STUBGEN_PACKAGE_ROOT)
        set(STUBGEN_PACKAGE_ROOT "")
    endif()

    # Locate Python
    set(Python3_ARTIFACTS_PREFIX "_HOST")
    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    # Run pybind11-stubgen in the installation prefix
    set(STUBGEN_MODULE ${STUBGEN_PACKAGE}.$<TARGET_FILE_BASE_NAME:${target}>)
    set(STUBGEN_CMD "\"${Python3_HOST_EXECUTABLE}\" -m pybind11_stubgen -o \"${ALPAQA_INSTALL_PYSTUBSDIR}\" 
        --exit-code --numpy-array-use-type-var --enum-class-locations Sign:LBFGS
        \"${STUBGEN_MODULE}\"")
    install(CODE "
        message(STATUS \"Executing pybind11-stubgen for ${STUBGEN_MODULE} \"
                       \"(destination: \\\"\${CMAKE_INSTALL_PREFIX}/${STUBGEN_PACKAGE_ROOT}\\\", interpreter: \\\"${Python3_HOST_EXECUTABLE}\\\")\")
        execute_process(COMMAND ${STUBGEN_CMD}
                        WORKING_DIRECTORY \"\${CMAKE_INSTALL_PREFIX}/${STUBGEN_PACKAGE_ROOT}\"
                        RESULT_VARIABLE STUBGEN_RET)
        if(NOT STUBGEN_RET EQUAL 0)
            message(SEND_ERROR \"pybind11-stubgen ${STUBGEN_MODULE} failed.\")
        endif()
        " EXCLUDE_FROM_ALL COMPONENT ${STUBGEN_COMPONENT})

endfunction()
