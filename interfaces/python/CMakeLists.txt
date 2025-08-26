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

# Find the Python development files
if (CMAKE_CROSSCOMPILING AND NOT DEFINED CMAKE_CROSSCOMPILING_EMULATOR)
    find_package(Python3 REQUIRED COMPONENTS Development.Module)
else()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
endif()

# Find the pybind11 library
set(PYBIND11_USE_CROSSCOMPILING On) # https://github.com/pybind/pybind11/pull/5083
find_package(pybind11 REQUIRED)

# Compile the example Python module
pybind11_add_module(_alpaqa MODULE
    "src/alpaqa.py.cpp"
    "src/accel/anderson.py.cpp"
    "src/accel/lbfgs.py.cpp"
    "src/outer/alm.py.cpp"
    "src/counters.py.cpp"
    "src/misc.py.cpp"
    "src/enums.py.cpp"
    "src/inner/panoc.py.cpp"
    "src/inner/fista.py.cpp"
    "src/inner/zerofpr.py.cpp"
    "src/inner/pantr.py.cpp"
    "src/inner/pantr-directions.py.cpp"
    "src/inner/panoc-directions.py.cpp"
    "src/inner/inner-solver.py.cpp"
    "src/inner/python-inner-solver.py.cpp"
    "src/problem/problems.py.cpp"
    "src/proximal/prox.py.cpp"
    # NO_EXTRAS # Prevent pybind11 from stripping the binary
)
if (ALPAQA_WITH_OCP)
    target_sources(_alpaqa PRIVATE
        "src/inner/panoc-ocp.py.cpp"
        "src/inner/ocp.py.cpp"
        "src/problem/control-problems.py.cpp")
endif()
if (ALPAQA_WITH_LBFGSB)
    target_sources(_alpaqa PRIVATE "src/inner/lbfgsb.py.cpp")
    target_link_libraries(_alpaqa PRIVATE alpaqa::lbfgsb-adapter)
endif()
if (ALPAQA_WITH_IPOPT)
    target_sources(_alpaqa PRIVATE "src/ipopt/ipopt.py.cpp")
    target_link_libraries(_alpaqa PRIVATE alpaqa::ipopt-adapter)
endif()
target_include_directories(_alpaqa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(_alpaqa PRIVATE pybind11::pybind11 alpaqa::alpaqa)
target_link_libraries(_alpaqa PRIVATE alpaqa::warnings)
target_compile_definitions(_alpaqa PRIVATE
    MODULE_NAME=$<TARGET_FILE_BASE_NAME:_alpaqa>)
set_target_properties(_alpaqa PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN true
    RELEASE_POSTFIX ""
    DEBUG_POSTFIX ""
    RELWITHDEBINFO_POSTFIX ""
    MINSIZEREL_POSTFIX ""
)
include(GenerateExportHeader)
generate_export_header(_alpaqa
    BASE_NAME alpaqa_python
    EXPORT_FILE_NAME export-python/alpaqa-python/export.h)
target_include_directories(_alpaqa PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export-python>)
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
    target_link_libraries(_alpaqa PRIVATE alpaqa::casadi-loader)
endif()
if (TARGET alpaqa::cutest-interface)
    target_link_libraries(_alpaqa PRIVATE alpaqa::cutest-interface)
endif()
if (TARGET alpaqa::dl-loader)
    target_link_libraries(_alpaqa PRIVATE alpaqa::dl-loader)
endif()
if (TARGET alpaqa::casadi-ocp-loader)
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
add_custom_target(python_modules DEPENDS _alpaqa)

if (DEFINED PY_BUILD_CMAKE_VERSION)
    set(ALPAQA_MODULE_INSTALL_DIR "")
else()
    set(ALPAQA_MODULE_INSTALL_DIR "${ALPAQA_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
endif()

# Install the Python source code
install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/${PY_BUILD_CMAKE_MODULE_NAME}
        EXCLUDE_FROM_ALL
        COMPONENT python_source
        DESTINATION "${ALPAQA_MODULE_INSTALL_DIR}"
        PATTERN "__pycache__" EXCLUDE)

# Install the Python module
install(TARGETS _alpaqa
        EXCLUDE_FROM_ALL
        COMPONENT python_modules
        DESTINATION ./${ALPAQA_MODULE_INSTALL_DIR}/${PY_BUILD_CMAKE_MODULE_NAME})

# Strip and install debug information
include(${PROJECT_SOURCE_DIR}/cmake/Debug.cmake)
alpaqa_install_debug_syms(_alpaqa python_modules_debug
    ./${ALPAQA_MODULE_INSTALL_DIR}/${PY_BUILD_CMAKE_MODULE_NAME}
    ./${ALPAQA_MODULE_INSTALL_DIR}/${PY_BUILD_CMAKE_MODULE_NAME})

# Generate stubs for the Python module
set(ALPAQA_WITH_PY_STUBS_DEFAULT On)
if (CMAKE_CROSSCOMPILING)
    set(ALPAQA_WITH_PY_STUBS_DEFAULT Off)
endif()
option(ALPAQA_WITH_PY_STUBS
    "Generate Python stub files (.pyi) for the Python module."
    ${ALPAQA_WITH_PY_STUBS_DEFAULT})
if (ALPAQA_WITH_PY_STUBS)
    include(cmake/Pybind11Stubgen.cmake)
    pybind11_stubgen(_alpaqa
        COMPONENT python_stubs
        PACKAGE ${PY_BUILD_CMAKE_MODULE_NAME}
        PACKAGE_ROOT "${ALPAQA_MODULE_INSTALL_DIR}")
endif()
