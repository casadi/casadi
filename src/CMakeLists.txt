# Target options
# --------------

add_library(panoc-alm SHARED
    "src/panoc.cpp"
    "src/box.cpp"
    "src/lbfgs.cpp"
    "src/alm.cpp"
)
add_library(panoc-alm::panoc-alm ALIAS panoc-alm)
set_target_properties(panoc-alm PROPERTIES EXPORT_NAME panoc-alm::panoc-alm)
target_include_directories(panoc-alm
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
target_link_libraries(panoc-alm PUBLIC Eigen3::Eigen)
# target_compile_definitions(panoc-alm PUBLIC PRINT_DEBUG_COUT=1)

# Interop
# -------

if (TARGET CUTEst::headers)
    add_library(cutest-loader STATIC src/interop/cutest/CUTEstLoader.cpp)
    target_link_libraries(cutest-loader 
        PRIVATE 
            CUTEst::headers
            ${CMAKE_DL_LIBS}
        PUBLIC 
            panoc-alm::panoc-alm)
    add_library(panoc-alm::cutest-loader ALIAS cutest-loader)
endif()

# Target installation
# -------------------

include(GNUInstallDirs)

set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/PANOC-ALM")

# Add the panoc-alm library to the "export-set", install the library files
install(TARGETS panoc-alm
    EXPORT PANOC-ALMTargets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" 
        COMPONENT lib
)

# Install the header files
INSTALL (DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT dev
    FILES_MATCHING PATTERN "*.h*"
)

# Install the export set for use with the install-tree
install(EXPORT PANOC-ALMTargets 
    FILE PANOC-ALMTargets.cmake
    DESTINATION "${INSTALL_CMAKE_DIR}" 
        COMPONENT dev
)

# Generate the config file that includes the exports
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in"
    "${PROJECT_BINARY_DIR}/PANOC-ALMConfig.cmake"
    INSTALL_DESTINATION "${INSTALL_CMAKE_DIR}"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/PANOC-ALMConfigVersion.cmake"
    VERSION "${PANOC-ALM_VERSION}"
    COMPATIBILITY SameMajorVersion
)

# Install the PANOC-ALMConfig.cmake and PANOC-ALMConfigVersion.cmake
install(FILES
    "${PROJECT_BINARY_DIR}/PANOC-ALMConfig.cmake"
    "${PROJECT_BINARY_DIR}/PANOC-ALMConfigVersion.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}" 
        COMPONENT dev)

# Add all targets to the build tree export set
export(EXPORT PANOC-ALMTargets
    FILE "${PROJECT_BINARY_DIR}/PANOC-ALMTargets.cmake")

# Export the package for use from the build tree:
# This registers the build tree with a global CMake-registry, so you can use
# find_package(PANOC-ALM) to find the package in the build tree
export(PACKAGE PANOC-ALM)
