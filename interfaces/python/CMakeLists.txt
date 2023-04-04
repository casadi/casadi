set(PY_FULL_VERSION ${PROJECT_VERSION}${PY_VERSION_SUFFIX})
set(PY_BUILD_CMAKE_MODULE_NAME "alpaqa" CACHE STRING "Python module name")
set(PY_BUILD_CMAKE_PACKAGE_NAME "alpaqa" CACHE STRING "Python package name")

# Make sure that the Python and CMake versions match
if (DEFINED PY_BUILD_CMAKE_PACKAGE_VERSION)
    if (NOT "${PY_BUILD_CMAKE_PACKAGE_VERSION}" MATCHES "^${PY_FULL_VERSION}$")
        message(FATAL_ERROR "Version number does not match "
                             "(${PY_BUILD_CMAKE_PACKAGE_VERSION} - ${PY_FULL_VERSION}).")
    endif()
endif()

# Find the Pybind11 headers
include(cmake/QueryPythonForPybind11.cmake)
find_pybind11_python_first()

# Compile the example Python module
pybind11_add_module(_alpaqa MODULE
    "src/alpaqa.py.cpp"
    "src/accel/anderson.py.cpp"
    "src/accel/lbfgs.py.cpp"
    "src/outer/alm.py.cpp"
    "src/counters.py.cpp"
    "src/enums.py.cpp"
    "src/inner/panoc.py.cpp"
    "src/inner/pantr.py.cpp"
    "src/inner/pantr-directions.py.cpp" # TODO
    "src/inner/panoc-directions.py.cpp"
    "src/inner/inner-solver.py.cpp"
    "src/problem/problems.py.cpp"
    "src/problem/control-problems.py.cpp"
    "src/params/alm-params.cpp"
    "src/params/anderson-direction-params.cpp"
    "src/params/anderson-params.cpp"
    "src/params/inner-solve-options.cpp"
    "src/params/lbfgs-direction-params.cpp"
    "src/params/lbfgs-params.cpp"
    "src/params/panoc-params.cpp"
    "src/params/structured-lbfgs-direction-params.cpp"
    "src/params/structured-newton-direction-params.cpp"
    "src/params/newton-tr-direction-params.cpp"
    "src/params/pantr-params.cpp"
    "src/params/steihaug-params.cpp"
    # NO_EXTRAS # Prevent pybind11 from stripping the binary
)
if (ALPAQA_WITH_OCP)
    target_sources(_alpaqa PRIVATE
        "src/inner/panoc-ocp.py.cpp"
        "src/inner/ocp.py.cpp"
        "src/params/panoc-ocp-params.cpp"
    )
endif()
target_include_directories(_alpaqa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(_alpaqa PRIVATE pybind11::pybind11 alpaqa::alpaqa)
target_link_libraries(_alpaqa PRIVATE alpaqa::warnings)
target_compile_definitions(_alpaqa PRIVATE
    MODULE_NAME=$<TARGET_FILE_BASE_NAME:_alpaqa>
    VERSION_INFO="${PY_FULL_VERSION}")
set_target_properties(_alpaqa PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN true
    RELEASE_POSTFIX ""
    DEBUG_POSTFIX ""
    RELWITHDEBINFO_POSTFIX ""
    MINSIZEREL_POSTFIX ""
)
set(ALPAQA_PYTHON_DEBUG_CONFIG "Debug" CACHE STRING "")
set_property(CACHE ALPAQA_PYTHON_DEBUG_CONFIG PROPERTY STRINGS
    ${CMAKE_CONFIGURATION_TYPES})
if (ALPAQA_PYTHON_DEBUG_CONFIG)
    string(TOUPPER ${ALPAQA_PYTHON_DEBUG_CONFIG} ALPAQA_PYTHON_DEBUG_CONFIG_UP)
    set_target_properties(_alpaqa PROPERTIES
        ${ALPAQA_PYTHON_DEBUG_CONFIG_UP}_POSTFIX "_d"
        PDB_NAME_${ALPAQA_PYTHON_DEBUG_CONFIG_UP} "_alpaqa_d")
endif()
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_options(_alpaqa PRIVATE "LINKER:--exclude-libs,ALL")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(_alpaqa PRIVATE "/bigobj")
endif()
if (TARGET alpaqa::casadi-loader)
    target_compile_definitions(_alpaqa PRIVATE ALPAQA_HAVE_CASADI)
    target_link_libraries(_alpaqa PRIVATE alpaqa::casadi-loader)
endif()
if (TARGET alpaqa::casadi-ocp-loader)
    target_compile_definitions(_alpaqa PRIVATE ALPAQA_HAVE_CASADI_OCP)
    target_link_libraries(_alpaqa PRIVATE alpaqa::casadi-ocp-loader)
endif()
if (ALPAQA_ENABLE_PCH)
    target_precompile_headers(_alpaqa PRIVATE
        <pybind11/chrono.h>
        <pybind11/eigen.h>
        <pybind11/functional.h>
        <pybind11/gil.h>
        <pybind11/iostream.h>
        <pybind11/stl.h>)
endif()

# Install the Python module
install(TARGETS _alpaqa
        EXCLUDE_FROM_ALL
        COMPONENT python_modules
        DESTINATION ${PY_BUILD_CMAKE_MODULE_NAME})

# Strip and install debug information
include(${PROJECT_SOURCE_DIR}/cmake/Debug.cmake)
alpaqa_install_debug_syms(_alpaqa python_modules_debug
    ${PY_BUILD_CMAKE_MODULE_NAME} ${PY_BUILD_CMAKE_MODULE_NAME})

# Generate stubs for the Python module
if (NOT CMAKE_CROSSCOMPILING)
    set(ALPAQA_WITH_PY_STUBS_DEFAULT On)
else()
    set(ALPAQA_WITH_PY_STUBS_DEFAULT Off)
endif()
option(ALPAQA_WITH_PY_STUBS
    "Generate Python stub files (.pyi) for the Python module."
    ${ALPAQA_WITH_PY_STUBS_DEFAULT})
if (ALPAQA_WITH_PY_STUBS)
    include(cmake/Pybind11Stubgen.cmake)
    pybind11_stubgen(_alpaqa)
    pybind11_stubgen_install(_alpaqa ${PY_BUILD_CMAKE_MODULE_NAME})
endif()
