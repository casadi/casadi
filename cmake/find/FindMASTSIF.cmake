find_package(SIFDecode REQUIRED)

# MASTSIF problems

find_path(MASTSIF_DIR NAMES mastsif.html sif.html
    HINTS
        ${MASTSIF}
        $ENV{MASTSIF}
        ${CUTEST}/../MASTSIF
        ${CUTEST}/../mastsif
        ${CUTEST}/../SIF
        ${CUTEST}/../sif
    CMAKE_FIND_ROOT_PATH_BOTH
)

set(MASTSIF_DIR ${MASTSIF_DIR} CACHE PATH "")
mark_as_advanced(MASTSIF_DIR)
find_package_handle_standard_args(MASTSIF
    REQUIRED_VARS
        MASTSIF_DIR
)

function(cutest_sif_problem PROBLEM_NAME)
    set(options "ALL")
    set(oneValueArgs SUFFIX)
    set(multiValueArgs OPTIONS)
    cmake_parse_arguments(PARSE_ARGV 1 CUTEST_SIF_PROBLEM
        "${options}" "${oneValueArgs}" "${multiValueArgs}")
    set(FULL_PROBLEM_NAME "${PROBLEM_NAME}${CUTEST_SIF_PROBLEM_SUFFIX}")
    if (NOT CUTEST_SIF_PROBLEM_ALL)
        set(CUTEST_SIF_PROBLEM_EXCLUDE_FROM_ALL EXCLUDE_FROM_ALL)
    endif()
    if (NOT TARGET CUTEst::problem-${FULL_PROBLEM_NAME})
        set(PROBLEM_DIR ${CMAKE_BINARY_DIR}/CUTEst/${FULL_PROBLEM_NAME})
        file(MAKE_DIRECTORY ${PROBLEM_DIR})
        add_custom_command(
            OUTPUT
                ${PROBLEM_DIR}/OUTSDIF.d
                ${PROBLEM_DIR}/AUTOMAT.d
                ${PROBLEM_DIR}/ELFUN.f
                ${PROBLEM_DIR}/EXTER.f
                ${PROBLEM_DIR}/GROUP.f
                ${PROBLEM_DIR}/RANGE.f
            COMMAND ${CMAKE_COMMAND} -E env 
                ARCHDEFS="${ARCHDEFS_DIR}" 
                SIFDECODE="${SIFDECODE_DIR}"
                MASTSIF="${MASTSIF_DIR}"
                MYARCH="${CUTEST_MYARCH}"
                "${SIFDECODE_EXE}"
                ${CUTEST_SIF_PROBLEM_OPTIONS}
                ${PROBLEM_NAME}
            MAIN_DEPENDENCY
                "${MASTSIF_DIR}/${PROBLEM_NAME}.SIF"
            WORKING_DIRECTORY
                ${PROBLEM_DIR}
            USES_TERMINAL
        )
        add_library(cutest-problem-${FULL_PROBLEM_NAME} MODULE
            ${CUTEST_SIF_PROBLEM_EXCLUDE_FROM_ALL}
            ${PROBLEM_DIR}/ELFUN.f
            ${PROBLEM_DIR}/EXTER.f
            ${PROBLEM_DIR}/GROUP.f
            ${PROBLEM_DIR}/RANGE.f
        )
        target_link_libraries(cutest-problem-${FULL_PROBLEM_NAME}
            PRIVATE
                CUTEst::cutest)
        set_target_properties(cutest-problem-${FULL_PROBLEM_NAME}
            PROPERTIES
                OUTPUT_NAME ${FULL_PROBLEM_NAME}
                LIBRARY_OUTPUT_DIRECTORY ${PROBLEM_DIR}
                PREFIX ""
                RELEASE_POSTFIX ""
                DEBUG_POSTFIX ""
                RELWITHDEBINFO_POSTFIX ""
                MINSIZEREL_POSTFIX "")
        add_library(CUTEst::problem-${FULL_PROBLEM_NAME} 
                    ALIAS cutest-problem-${FULL_PROBLEM_NAME})
        add_custom_command(TARGET cutest-problem-${FULL_PROBLEM_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND}
                ARGS -E copy "${PROBLEM_DIR}/OUTSDIF.d"
                    "$<TARGET_FILE_DIR:cutest-problem-${FULL_PROBLEM_NAME}>")
    endif()
endfunction()
