include(GNUInstallDirs)
include(${PROJECT_SOURCE_DIR}/cmake/Debug.cmake)

set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/alpaqa")

# Add the alpaqa library to the "export-set", install the library files
install(TARGETS ${ALPAQA_INSTALL_TARGETS} warnings
    EXPORT alpaqaTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        COMPONENT shlib
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        COMPONENT shlib
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" 
        COMPONENT lib)

# Install the header files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/alpaqa/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT dev
    FILES_MATCHING REGEX "/.*\.(h|[hti]pp)$")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/interop/casadi/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT dev
    FILES_MATCHING REGEX "/.*\.(h|[hti]pp)$")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/interop/dl/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT dev
    FILES_MATCHING REGEX "/.*\.(h|[hti]pp)$")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/interop/dl-api/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT dl_dev
    FILES_MATCHING REGEX "/.*\.(h|[hti]pp)$")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/export/alpaqa/export.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/alpaqa/"
        COMPONENT dev)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/export/alpaqa/casadi-loader-export.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/alpaqa/"
        COMPONENT dev)

# Install the debug files
foreach(TGT IN LISTS ALPAQA_INSTALL_TARGETS)
    get_target_property(TGT_TYPE ${TGT} TYPE)
    if (${TGT_TYPE} STREQUAL "SHARED_LIBRARY")
        alpaqa_install_debug_syms(${TGT} debug
                                  ${CMAKE_INSTALL_LIBDIR}
                                  ${CMAKE_INSTALL_BINDIR})
    endif()
endforeach()

# Install the export set for use with the install tree
install(EXPORT alpaqaTargets 
    FILE alpaqaTargets.cmake
    DESTINATION "${INSTALL_CMAKE_DIR}" 
        COMPONENT dev
    NAMESPACE alpaqa::)

# Generate the config file that includes the exports
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in"
    "${PROJECT_BINARY_DIR}/alpaqaConfig.cmake"
    INSTALL_DESTINATION "${INSTALL_CMAKE_DIR}"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/alpaqaConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion)

# Install the alpaqaConfig.cmake and alpaqaConfigVersion.cmake
install(FILES
    "${PROJECT_BINARY_DIR}/alpaqaConfig.cmake"
    "${PROJECT_BINARY_DIR}/alpaqaConfigVersion.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}" 
        COMPONENT dev)

# Add all targets to the build tree export set
export(EXPORT alpaqaTargets
    FILE "${PROJECT_BINARY_DIR}/alpaqaTargets.cmake"
    NAMESPACE alpaqa::)
