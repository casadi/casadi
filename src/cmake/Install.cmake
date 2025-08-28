include(${PROJECT_SOURCE_DIR}/cmake/Debug.cmake)

# Set the runtime linker/loader search paths to make alpaqa stand-alone
if (ALPAQA_STANDALONE)
    cmake_path(RELATIVE_PATH ALPAQA_INSTALL_LIBDIR
               BASE_DIRECTORY ALPAQA_INSTALL_BINDIR
               OUTPUT_VARIABLE ALPAQA_INSTALL_LIBRELBINDIR)
    foreach (TGT IN LISTS ALPAQA_INSTALL_TARGETS ALPAQA_INSTALL_EXE)
        if (APPLE)
            set_target_properties(${TGT} PROPERTIES
                INSTALL_RPATH "$@loader_path;$@loader_path/${ALPAQA_INSTALL_LIBRELBINDIR}")
        else()
            set_target_properties(${TGT} PROPERTIES
                INSTALL_RPATH "$ORIGIN;$ORIGIN/${ALPAQA_INSTALL_LIBRELBINDIR}")
        endif()
    endforeach()
endif()

function(alpaqa_add_if_target_exists OUT)
    foreach(TGT IN LISTS ARGN)
        if (TARGET ${TGT})
            list(APPEND ${OUT} ${TGT})
        endif()
    endforeach()
    set(${OUT} ${${OUT}} PARENT_SCOPE)
endfunction()

include(CMakePackageConfigHelpers)

set(ALPAQA_INSTALLED_COMPONENTS)
macro(alpaqa_install_config PKG COMP)
    # Install the target CMake definitions
    install(EXPORT alpaqa${PKG}Targets
        FILE alpaqa${PKG}Targets.cmake
        DESTINATION "${ALPAQA_INSTALL_CMAKEDIR}"
            COMPONENT ${COMP}
        NAMESPACE alpaqa::)
    # Add all targets to the build tree export set
    export(EXPORT alpaqa${PKG}Targets
        FILE "${PROJECT_BINARY_DIR}/alpaqa${PKG}Targets.cmake"
        NAMESPACE alpaqa::)
    # Generate the config file that includes the exports
    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PKG}Config.cmake.in"
        "${PROJECT_BINARY_DIR}/alpaqa${PKG}Config.cmake"
        INSTALL_DESTINATION "${ALPAQA_INSTALL_CMAKEDIR}"
        NO_SET_AND_CHECK_MACRO)
    write_basic_package_version_file(
        "${PROJECT_BINARY_DIR}/alpaqa${PKG}ConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}"
        COMPATIBILITY SameMinorVersion)
    # Install the alpaqaConfig.cmake and alpaqaConfigVersion.cmake
    install(FILES
        "${PROJECT_BINARY_DIR}/alpaqa${PKG}Config.cmake"
        "${PROJECT_BINARY_DIR}/alpaqa${PKG}ConfigVersion.cmake"
        DESTINATION "${ALPAQA_INSTALL_CMAKEDIR}"
            COMPONENT ${COMP})
    list(APPEND ALPAQA_OPTIONAL_COMPONENTS ${PKG})
endmacro()

macro(alpaqa_install_headers DIR COMP)
    # Install the header files
    install(DIRECTORY ${DIR}
        DESTINATION "${ALPAQA_INSTALL_INCLUDEDIR}"
            COMPONENT ${COMP}
        FILES_MATCHING REGEX "/.*\\.(h|[hti]pp)$")
endmacro()

macro(alpaqa_install_cmake FILES COMP)
    # Install a CMake script
    install(FILES ${FILES}
        DESTINATION "${ALPAQA_INSTALL_CMAKEDIR}"
            COMPONENT ${COMP})
endmacro()

set(ALPAQA_INSTALLED_TARGETS_MSG "\nSummary of alpaqa components and targets to install:\n\n")

# Install the alpaqa core libraries
install(TARGETS warnings alpaqa
    EXPORT alpaqaCoreTargets
    RUNTIME DESTINATION "${ALPAQA_INSTALL_BINDIR}"
        COMPONENT lib
    LIBRARY DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
        COMPONENT lib
        NAMELINK_COMPONENT dev
    ARCHIVE DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
        COMPONENT dev)
alpaqa_install_config(Core dev)
alpaqa_install_headers("${PROJECT_BINARY_DIR}/include/" dev)
alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/alpaqa/include/" dev)
alpaqa_install_headers("${CMAKE_CURRENT_BINARY_DIR}/export/" dev)
string(APPEND ALPAQA_INSTALLED_TARGETS_MSG " * Core:   alpaqa\n")

# Install the CasADi interface
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_CASADI_TARGETS "casadi-loader")
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_CASADI_TARGETS "casadi-ocp-loader")
if (ALPAQA_COMPONENT_CASADI_TARGETS)
    install(TARGETS ${ALPAQA_COMPONENT_CASADI_TARGETS}
        EXPORT alpaqaCasADiTargets
        RUNTIME DESTINATION "${ALPAQA_INSTALL_BINDIR}"
            COMPONENT casadi
        LIBRARY DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
            COMPONENT casadi
            NAMELINK_COMPONENT casadi_dev
        ARCHIVE DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
            COMPONENT casadi_dev)
    alpaqa_install_config(CasADi casadi_dev)
    alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/interop/casadi/include/" casadi_dev)
    list(JOIN ALPAQA_COMPONENT_CASADI_TARGETS ", " TGTS)
    string(APPEND ALPAQA_INSTALLED_TARGETS_MSG " * CasADi: ${TGTS}\n")
endif()

# Install the DL API
install(TARGETS dl-api
    EXPORT alpaqaDlTargets)
alpaqa_install_config(Dl dl_dev)
alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/interop/dl-api/include/" dl_dev)
alpaqa_install_cmake("${CMAKE_CURRENT_SOURCE_DIR}/cmake/dl-problem.cmake" dl_dev)
string(APPEND ALPAQA_INSTALLED_TARGETS_MSG " * Dl:     dl-api\n")

# Install everything else
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_EXTRA_TARGETS "dl-loader")
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_EXTRA_TARGETS "cutest-interface")
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_EXTRA_TARGETS "ipopt-adapter")
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_EXTRA_TARGETS "lbfgsb-fortran")
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_EXTRA_TARGETS "lbfgsb-adapter")
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_EXTRA_TARGETS "qpalm-adapter")
if (ALPAQA_COMPONENT_EXTRA_TARGETS)
    install(TARGETS ${ALPAQA_COMPONENT_EXTRA_TARGETS}
        EXPORT alpaqaExtraTargets
        RUNTIME DESTINATION "${ALPAQA_INSTALL_BINDIR}"
            COMPONENT extra
        LIBRARY DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
            COMPONENT extra
            NAMELINK_COMPONENT extra_dev
        ARCHIVE DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
            COMPONENT extra_dev)
    alpaqa_install_config(Extra extra_dev)
    alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/interop/cutest/include/" extra_dev)
    alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/interop/ipopt/include/" extra_dev)
    alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/interop/lbfgsb/include/" extra_dev)
    alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/interop/dl/include/" extra_dev)
    alpaqa_install_headers("${PROJECT_SOURCE_DIR}/src/interop/qpalm/include/" extra_dev)
    list(JOIN ALPAQA_COMPONENT_EXTRA_TARGETS ", " TGTS)
    string(APPEND ALPAQA_INSTALLED_TARGETS_MSG " * Extra:  ${TGTS}\n")
endif()

# Install the tools
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_TOOLS_TARGETS "driver")
alpaqa_add_if_target_exists(ALPAQA_COMPONENT_TOOLS_TARGETS "gradient-checker")
if (ALPAQA_COMPONENT_TOOLS_TARGETS)
    install(TARGETS ${ALPAQA_COMPONENT_TOOLS_TARGETS}
        EXPORT alpaqaToolsTargets
        RUNTIME DESTINATION "${ALPAQA_INSTALL_BINDIR}"
            COMPONENT bin
        LIBRARY DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
            COMPONENT bin
            NAMELINK_COMPONENT bin
        ARCHIVE DESTINATION "${ALPAQA_INSTALL_LIBDIR}"
            COMPONENT bin)
    alpaqa_install_config(Tools bin)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/alpaqa/src/driver/alpaqa-driver.autocomplete.zsh"
        DESTINATION "${ALPAQA_INSTALL_ZSHCOMPLETEDIR}/"
        PERMISSIONS OWNER_READ;OWNER_WRITE;OWNER_EXECUTE;GROUP_READ;GROUP_EXECUTE;WORLD_READ;WORLD_EXECUTE
        COMPONENT bin)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/alpaqa/src/driver/alpaqa-driver.autocomplete.bash"
        DESTINATION "${ALPAQA_INSTALL_BASHCOMPLETEDIR}/"
        RENAME "alpaqa-driver"
        PERMISSIONS OWNER_READ;OWNER_WRITE;GROUP_READ;WORLD_READ
        COMPONENT bin)
    list(JOIN ALPAQA_COMPONENT_TOOLS_TARGETS ", " TGTS)
    string(APPEND ALPAQA_INSTALLED_TARGETS_MSG " * Tools:  ${TGTS}\n")
endif()

# Install the debug files
foreach(target IN LISTS ALPAQA_INSTALL_TARGETS ALPAQA_INSTALL_EXE)
    get_target_property(target_type ${target} TYPE)
    if (${target_type} STREQUAL "SHARED_LIBRARY")
        alpaqa_install_debug_syms(${target} debug
                                  ${ALPAQA_INSTALL_LIBDIR}
                                  ${ALPAQA_INSTALL_BINDIR})
    elseif (${target_type} STREQUAL "EXECUTABLE")
        alpaqa_install_debug_syms(${target} debug
                                  ${ALPAQA_INSTALL_BINDIR}
                                  ${ALPAQA_INSTALL_BINDIR})
    endif()
endforeach()

# Generate the main config file
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in"
    "${PROJECT_BINARY_DIR}/alpaqaConfig.cmake"
    INSTALL_DESTINATION "${ALPAQA_INSTALL_CMAKEDIR}"
    NO_SET_AND_CHECK_MACRO)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/alpaqaConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMinorVersion)
# Install the main alpaqaConfig.cmake and alpaqaConfigVersion.cmake files
install(FILES
    "${PROJECT_BINARY_DIR}/alpaqaConfig.cmake"
    "${PROJECT_BINARY_DIR}/alpaqaConfigVersion.cmake"
    DESTINATION "${ALPAQA_INSTALL_CMAKEDIR}"
        COMPONENT dev)

# Print the components and targets we're going to install
message(${ALPAQA_INSTALLED_TARGETS_MSG})