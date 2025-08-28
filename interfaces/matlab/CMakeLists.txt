set(ALPAQA_INSTALL_MEXDIR "${ALPAQA_INSTALL_LIBDIR}/alpaqa/mex" CACHE PATH
    "Installation directory for the alpaqa MATLAB/MEX interface")

# Matlab uses UTF-16 for all its strings, so we need to convert back and forth
find_package(utf8cpp REQUIRED)
if (NOT TARGET utf8cpp::utf8cpp) # Conan target is named differently
    add_library(utf8cpp::utf8cpp ALIAS utf8cpp)
endif()

# Matlab interface requires CasADi and JSON
if (NOT TARGET alpaqa::casadi-loader)
    message(FATAL_ERROR "MATLAB MEX interface requires CasADi")
endif()
if (NOT ALPAQA_WITH_JSON)
    message(FATAL_ERROR "MATLAB MEX interface requires JSON")
endif()

# Find the Matlab MEX compiler
find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER)

# Create a target for compiling a MEX function
matlab_add_mex(NAME alpaqa_mex
    SRC
        "src/alpaqa_mex.cpp" "src/minimize.cpp" "src/unicode.cpp"
    LINK_TO
        alpaqa::alpaqa alpaqa::casadi-loader utf8cpp::utf8cpp alpaqa::warnings
)
target_include_directories(alpaqa_mex PRIVATE include)
alpaqa_configure_visibility(alpaqa_mex)
if (NOT WIN32 AND NOT APPLE)
    # Matlab's C++ standard library is too old, so link it statically
    target_link_options(alpaqa_mex PUBLIC -static-libstdc++)
endif()
# Get MSVC to report the correct value for __cplusplus
if (MSVC)
    target_compile_options(alpaqa_mex PRIVATE /Zc:__cplusplus)
endif()

# Don't use any filename postfixes for the MEX file
set_target_properties(alpaqa_mex PROPERTIES
        RELEASE_POSTFIX ""
        DEBUG_POSTFIX ""
        RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX "")

# Set the runtime linker/loader search paths to make alpaqa stand-alone
if (ALPAQA_STANDALONE)
    cmake_path(RELATIVE_PATH ALPAQA_INSTALL_LIBDIR
               BASE_DIRECTORY "${ALPAQA_INSTALL_MEXDIR}/+alpaqa"
               OUTPUT_VARIABLE ALPAQA_INSTALL_LIBRELMEXDIR)
    if (APPLE)
        set_target_properties(alpaqa_mex PROPERTIES INSTALL_RPATH
            "@loader_path;@loader_path/${ALPAQA_INSTALL_LIBRELMEXDIR}")
    else()
        set_target_properties(alpaqa_mex PROPERTIES INSTALL_RPATH
            "$ORIGIN;$ORIGIN/${ALPAQA_INSTALL_LIBRELMEXDIR}")
    endif()
endif()

# Install the final MEX file
set(ALPAQA_MATLAB_FILES
    "minimize.m"
    "version.m"
    "Problem.m"
)
# Install to root of install prefix
install(TARGETS alpaqa_mex
        EXCLUDE_FROM_ALL
        COMPONENT mex_interface
        DESTINATION "+alpaqa")
install(FILES ${ALPAQA_MATLAB_FILES}
        EXCLUDE_FROM_ALL
        COMPONENT mex_interface
        DESTINATION "+alpaqa")
# Install into lib/alpaqa/mex subdirectories (default)
install(TARGETS alpaqa_mex
        COMPONENT matlab
        DESTINATION "${ALPAQA_INSTALL_MEXDIR}/+alpaqa")
install(FILES ${ALPAQA_MATLAB_FILES}
        COMPONENT matlab
        DESTINATION "${ALPAQA_INSTALL_MEXDIR}/+alpaqa")
