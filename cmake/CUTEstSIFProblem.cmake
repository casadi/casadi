macro(cutest_sif_problem PROBLEM_NAME)
    set(PROBLEM_DIR CUTEst/${PROBLEM_NAME})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PROBLEM_DIR})
    add_custom_command(
        OUTPUT
            ${PROBLEM_DIR}/OUTSDIF.d
            ${PROBLEM_DIR}/AUTOMAT.d
            ${PROBLEM_DIR}/ELFUN.f
            ${PROBLEM_DIR}/EXTER.f
            ${PROBLEM_DIR}/GROUP.f
            ${PROBLEM_DIR}/RANGE.f
        COMMAND ${CMAKE_COMMAND} -E env 
            ARCHDEFS="${ARCHDefs_DIR}" 
            SIFDECODE="${SIFDecoder_DIR}"
            MASTSIF="${MASTSIF_DIR}"
            MYARCH="${CUTEst_MYARCH}"
            "${SIFDecoder_EXE}"
            ${PROBLEM_NAME}
        MAIN_DEPENDENCY
            "${MASTSIF_DIR}/${PROBLEM_NAME}.SIF"
        WORKING_DIRECTORY
            ${PROBLEM_DIR}
    )
    add_library(CUTEst_${PROBLEM_NAME} SHARED 
        ${PROBLEM_DIR}/ELFUN.f
        ${PROBLEM_DIR}/EXTER.f
        ${PROBLEM_DIR}/GROUP.f
        ${PROBLEM_DIR}/RANGE.f
    )
    target_link_libraries(CUTEst_${PROBLEM_NAME} PRIVATE CUTEst::objects)
    set_target_properties(CUTEst_${PROBLEM_NAME}
        PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${PROBLEM_DIR}
            DEBUG_POSTFIX ""
            ASAN_POSTFIX ""
            TSAN_POSTFIX "")
endmacro()
