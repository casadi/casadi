set(PY_FULL_VERSION ${PROJECT_VERSION}${PY_VERSION_SUFFIX})

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
Python3_add_library(_alpaqa MODULE "src/alpaqa.py.cpp")
target_link_libraries(_alpaqa PRIVATE pybind11::pybind11 alpaqa::alpaqa)
target_compile_definitions(_alpaqa PRIVATE
    MODULE_NAME=$<TARGET_FILE_BASE_NAME:_alpaqa>
    VERSION_INFO="${PY_FULL_VERSION}")
set_target_properties(_alpaqa PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN true)
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_options(_alpaqa PRIVATE "LINKER:--exclude-libs,ALL")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(_alpaqa PRIVATE "/bigobj")
endif()
if (TARGET alpaqa::casadi-loader)
    target_compile_definitions(_alpaqa PRIVATE ALPAQA_HAVE_CASADI)
    target_link_libraries(_alpaqa PRIVATE alpaqa::casadi-loader)
endif()

# Install the Python module
install(TARGETS _alpaqa
        EXCLUDE_FROM_ALL
        COMPONENT python_modules
        DESTINATION ./${PY_BUILD_CMAKE_MODULE_NAME})

# Generate stubs for the Python module
option(ALPAQA_WITH_PY_STUBS
    "Generate Python stub files (.pyi) for the Python module." On)
if (ALPAQA_WITH_PY_STUBS AND NOT CMAKE_CROSSCOMPILING)
    include(cmake/Pybind11Stubgen.cmake)
    pybind11_stubgen(_alpaqa)
    pybind11_stubgen_install(_alpaqa ./${PY_BUILD_CMAKE_MODULE_NAME})
endif()
