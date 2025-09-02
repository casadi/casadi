function(alpaqa_configure_dl_problem_visibility target)
    cmake_parse_arguments(ALPAQA_CONFIG_VIS "" "FUNCTION_NAME" "" ${ARGN})
    if (NOT DEFINED ALPAQA_CONFIG_VIS_FUNCTION_NAME)
        set(ALPAQA_CONFIG_VIS_FUNCTION_NAME "register_alpaqa_problem")
    endif()
    set_target_properties(${target} PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                               C_VISIBILITY_PRESET "hidden"
                                               VISIBILITY_INLINES_HIDDEN true)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        set(VERSION_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/${target}-export.lds")
        file(WRITE ${VERSION_SCRIPT}
            "{\n"
            "  global:\n"
            "    ${ALPAQA_CONFIG_VIS_FUNCTION_NAME};\n"
            "    ${ALPAQA_CONFIG_VIS_FUNCTION_NAME}_version;\n"
            "    _ZTI*;\n"
            "    _ZTS*;\n\n"
            "  local:\n"
            "    *;\n"
            "};")
        target_link_options(${target} PRIVATE
            "LINKER:--version-script=${VERSION_SCRIPT}"
            "LINKER:--exclude-libs,ALL")
    endif()
endfunction()

#[=======================================================================[.rst:

.. cmake:command::
    alpaqa_add_dl_problem_module

    Create a target for a problem that can be dymanically loaded by the alpaqa
    solvers.

    .. code-block:: cmake

        alpaqa_add_dl_problem_module(<target> [LINK_ALPAQA] [VISIBILITY_DEFAULT] [NO_CXX_20]
                                     [FUNCTION_NAME <name>] [FILES [<file> ...]])

    This function creates a CMake module library target with the given name
    ``<target>``, configures the correct visibility settings, generates an
    export header, and links to the correct alpaqa targets.
    
    If the ``FILES`` option is given, the provided files are added to the
    target. Otherwise, the source file ``${target}.cpp`` in the current folder
    is used.

    In addition to linking to the ``alpaqa::dl-api`` headers, the main
    ``alpaqa::alpaqa`` library is linked as well if the ``LINK_ALPAQA``
    option is provided. This is useful if you need any of the utility functions
    provided by alpaqa.

    If the ``VISIBILITY_DEFAULT`` option is given the visibility settings are
    not altered.

    By default, C++20 support is enabled, unless the ``NO_CXX_20`` is given.

    ``FUNCTION_NAME`` can be used to set the name of the registration function
    to expose. It should match the function defined in the module (which
    should have C linkage, using ``extern "C"``).

    .. note::

        This function requires the alpaqa ``Dl`` component to be loaded:

        .. code-block:: cmake

            find_package(alpaqa 1.1.0 REQUIRED COMPONENTS Dl)

        If ``LINK_ALPAQA`` is specified, both the ``Core`` and ``Dl`` components
        are required.

#]=======================================================================]
function(alpaqa_add_dl_problem_module target)
    cmake_parse_arguments(ALPAQA_PROBLEM_MODULE
        "LINK_ALPAQA;VISIBILITY_DEFAULT;NO_CXX_20"
        "FUNCTION_NAME"
        "FILES" ${ARGN})
    # The resulting binary is meant to be loaded dynamically, so we use the
    # MODULE type (not SHARED).
    if (DEFINED ALPAQA_PROBLEM_MODULE_FILES)
        add_library(${target} MODULE ${ALPAQA_PROBLEM_MODULE_FILES})
    else()
        add_library(${target} MODULE "${target}.cpp")
    endif()
    # We only depend on the dl-api header, but the user might also want to link
    # to the main alpaqa library for some convenience functions.
    target_link_libraries(${target} PRIVATE alpaqa::dl-api)
    if (ALPAQA_PROBLEM_MODULE_LINK_ALPAQA)
        if (NOT TARGET alpaqa::alpaqa)
            message(FATAL_ERROR "The LINK_ALPAQA option requires the alpaqa Core component to be loaded. "
                "For example, use find_package(alpaqa ${ALPAQA_VERSION} REQUIRED COMPONENTS Core Dl).\n")
        endif()
        target_link_libraries(${target} PRIVATE alpaqa::alpaqa)
    endif()
    # Enable C++20 for optional features of the DL api
    if (NOT ALPAQA_PROBLEM_MODULE_NO_CXX_20)
        target_compile_features(${target} PRIVATE cxx_std_20)
    endif()
    # By default, we want to hide all symbols except the main entry point, to
    # avoid conflicts and to keep the interface clean
    if (NOT ALPAQA_PROBLEM_MODULE_VISIBILITY_DEFAULT)
        if (DEFINED ALPAQA_PROBLEM_MODULE_FUNCTION_NAME)
            alpaqa_configure_dl_problem_visibility(${target}
                FUNCTION_NAME ${ALPAQA_PROBLEM_MODULE_FUNCTION_NAME})
        else()
            alpaqa_configure_dl_problem_visibility(${target})
        endif()
    endif()
    # We don't want any prefixes or postfixes to modify the file name
    set_target_properties(${target} PROPERTIES
        PREFIX "" RELEASE_POSTFIX "" DEBUG_POSTFIX "" RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX "")
    # Finally, we generate a header with macros to export the entry point
    # symbols (e.g. using [[gnu::visibility("default")]] for GCC and Clang, or
    # __declspec(dllexport) for MSVC)
    include(GenerateExportHeader)
    generate_export_header(${target}
        EXPORT_FILE_NAME export-${target}/${target}/export.h)
    target_include_directories(${target} PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export-${target}>)
endfunction()
