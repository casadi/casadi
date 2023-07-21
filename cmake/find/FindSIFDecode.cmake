find_package(CUTEst REQUIRED)

# SIFDecode

find_program(SIFDECODE_EXE sifdecoder
    HINTS
        ${SIFDECODE}
        $ENV{SIFDECODE}
        ${CUTEST}/../SIFDecode
        ${CUTEST}/../sifdecode
    PATH_SUFFIXES
        bin
    NO_CMAKE_FIND_ROOT_PATH
)
cmake_path(GET SIFDECODE_EXE PARENT_PATH SIFDECODE_EXE_DIR)
set(SIFDECODE_DIR ${SIFDECODE_EXE_DIR}/.. CACHE PATH "")

find_path(ARCHDEFS_DIR NAMES compiler.${CUTEST_MYARCH}
    HINTS
        ${ARCHDEFS}
        $ENV{ARCHDEFS}
        ${CUTEST}/../ARCHDefs
        ${CUTEST}/../archdefs
)
set(ARCHDEFS ${ARCHDEFS_DIR} CACHE PATH "")

mark_as_advanced(SIFDECODE_EXE ARCHDEFS_DIR)
find_package_handle_standard_args(SIFDecode 
    REQUIRED_VARS
        SIFDECODE_EXE
        ARCHDEFS_DIR
)
