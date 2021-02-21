set(CUTEst_MYARCH "pc64.lnx.gfo" CACHE STRING "CUTEst architecture" )

# CUTEst itself

find_path(CUTEst_DIR include/cutest.h
    HINTS
        ${CUTEST}
        ENV CUTEST
    PATH_SUFFIXES
        CUTEst
        cutest
        opt
        opt/CUTEst
        opt/cutest
        opt/CUTEst/cutest
        opt/cutest/cutest
)
find_library(CUTEst_LIBRARY
    NAMES 
        cutest
    PATHS
        ${CUTEst_DIR}
    PATH_SUFFIXES
        objects/${CUTEst_MYARCH}/double
)
mark_as_advanced(CUTEst_DIR CUTEst_LIBRARY)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUTEst
    REQUIRED_VARS
        CUTEst_DIR
        CUTEst_LIBRARY
)

# CUTEst library
if(CUTEst_FOUND AND NOT TARGET CUTEst::cutest)
    add_library(CUTEst::cutest STATIC IMPORTED)
    set_target_properties(CUTEst::cutest PROPERTIES
        IMPORTED_LOCATION ${CUTEst_LIBRARY})
    target_include_directories(CUTEst::cutest
        INTERFACE ${CUTEst_DIR}/include)
endif()
if(CUTEst_FOUND AND NOT TARGET CUTEst::objects)
    execute_process(COMMAND ${CMAKE_AR} x ${CUTEst_LIBRARY}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set(CUTEst_OBJECTS
        ${CMAKE_CURRENT_BINARY_DIR}/cutest.o
        ${CMAKE_CURRENT_BINARY_DIR}/pname.o
        ${CMAKE_CURRENT_BINARY_DIR}/probname.o
        ${CMAKE_CURRENT_BINARY_DIR}/varnames.o
        ${CMAKE_CURRENT_BINARY_DIR}/newthread.o
        ${CMAKE_CURRENT_BINARY_DIR}/problem.o
        ${CMAKE_CURRENT_BINARY_DIR}/fortran_ops.o
        ${CMAKE_CURRENT_BINARY_DIR}/interface.o
        ${CMAKE_CURRENT_BINARY_DIR}/ccutest.o
        ${CMAKE_CURRENT_BINARY_DIR}/timings.o
        ${CMAKE_CURRENT_BINARY_DIR}/usetup.o
        ${CMAKE_CURRENT_BINARY_DIR}/udimen.o
        ${CMAKE_CURRENT_BINARY_DIR}/udimse.o
        ${CMAKE_CURRENT_BINARY_DIR}/udimsh.o
        ${CMAKE_CURRENT_BINARY_DIR}/unames.o
        ${CMAKE_CURRENT_BINARY_DIR}/uvartype.o
        ${CMAKE_CURRENT_BINARY_DIR}/ufn.o
        ${CMAKE_CURRENT_BINARY_DIR}/ugr.o
        ${CMAKE_CURRENT_BINARY_DIR}/uofg.o
        ${CMAKE_CURRENT_BINARY_DIR}/udh.o
        ${CMAKE_CURRENT_BINARY_DIR}/ugrdh.o
        ${CMAKE_CURRENT_BINARY_DIR}/ush.o
        ${CMAKE_CURRENT_BINARY_DIR}/ushp.o
        ${CMAKE_CURRENT_BINARY_DIR}/ueh.o
        ${CMAKE_CURRENT_BINARY_DIR}/ugreh.o
        ${CMAKE_CURRENT_BINARY_DIR}/ugrsh.o
        ${CMAKE_CURRENT_BINARY_DIR}/uhprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/ushprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/ubandh.o
        ${CMAKE_CURRENT_BINARY_DIR}/ureport.o
        ${CMAKE_CURRENT_BINARY_DIR}/uterminate.o
        ${CMAKE_CURRENT_BINARY_DIR}/csetup.o
        ${CMAKE_CURRENT_BINARY_DIR}/cdimen.o
        ${CMAKE_CURRENT_BINARY_DIR}/cdimse.o
        ${CMAKE_CURRENT_BINARY_DIR}/cdimsh.o
        ${CMAKE_CURRENT_BINARY_DIR}/cdimsj.o
        ${CMAKE_CURRENT_BINARY_DIR}/cdimchp.o
        ${CMAKE_CURRENT_BINARY_DIR}/cnames.o
        ${CMAKE_CURRENT_BINARY_DIR}/cvartype.o
        ${CMAKE_CURRENT_BINARY_DIR}/cfn.o
        ${CMAKE_CURRENT_BINARY_DIR}/cgr.o
        ${CMAKE_CURRENT_BINARY_DIR}/cofg.o
        ${CMAKE_CURRENT_BINARY_DIR}/cofsg.o
        ${CMAKE_CURRENT_BINARY_DIR}/ccfg.o
        ${CMAKE_CURRENT_BINARY_DIR}/clfg.o
        ${CMAKE_CURRENT_BINARY_DIR}/ccfsg.o
        ${CMAKE_CURRENT_BINARY_DIR}/ccifg.o
        ${CMAKE_CURRENT_BINARY_DIR}/ccifsg.o
        ${CMAKE_CURRENT_BINARY_DIR}/cdh.o
        ${CMAKE_CURRENT_BINARY_DIR}/cdhc.o
        ${CMAKE_CURRENT_BINARY_DIR}/ceh.o
        ${CMAKE_CURRENT_BINARY_DIR}/cgrdh.o
        ${CMAKE_CURRENT_BINARY_DIR}/cifn.o
        ${CMAKE_CURRENT_BINARY_DIR}/cigr.o
        ${CMAKE_CURRENT_BINARY_DIR}/cisgr.o
        ${CMAKE_CURRENT_BINARY_DIR}/cidh.o
        ${CMAKE_CURRENT_BINARY_DIR}/csh.o
        ${CMAKE_CURRENT_BINARY_DIR}/cshc.o
        ${CMAKE_CURRENT_BINARY_DIR}/cshp.o
        ${CMAKE_CURRENT_BINARY_DIR}/cish.o
        ${CMAKE_CURRENT_BINARY_DIR}/cjprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/cstats.o
        ${CMAKE_CURRENT_BINARY_DIR}/csgr.o
        ${CMAKE_CURRENT_BINARY_DIR}/csgreh.o
        ${CMAKE_CURRENT_BINARY_DIR}/csgrsh.o
        ${CMAKE_CURRENT_BINARY_DIR}/csjprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/chprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/chcprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/cshprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/cshcprod.o
        ${CMAKE_CURRENT_BINARY_DIR}/cchprods.o
        ${CMAKE_CURRENT_BINARY_DIR}/csjp.o
        ${CMAKE_CURRENT_BINARY_DIR}/csgrp.o
        ${CMAKE_CURRENT_BINARY_DIR}/csgrshp.o
        ${CMAKE_CURRENT_BINARY_DIR}/cchprodsp.o
        ${CMAKE_CURRENT_BINARY_DIR}/creport.o
        ${CMAKE_CURRENT_BINARY_DIR}/connames.o
        ${CMAKE_CURRENT_BINARY_DIR}/cterminate.o
        ${CMAKE_CURRENT_BINARY_DIR}/lqp.o
        ${CMAKE_CURRENT_BINARY_DIR}/cconst.o)
    add_library(CUTEst::objects OBJECT IMPORTED)
    set_target_properties(CUTEst::objects PROPERTIES
        LINKER_LANGUAGE Fortran
        IMPORTED_OBJECTS "${CUTEst_OBJECTS}")
endif()
if(CUTEst_FOUND AND NOT TARGET CUTEst::headers)
    add_library(CUTEst::headers INTERFACE IMPORTED)
    target_include_directories(CUTEst::headers
        INTERFACE ${CUTEst_DIR}/include)
endif()

# SIFDecoder

find_path(SIFDecoder_DIR bin/sifdecoder
    HINTS
        ${SIFDECODE}
        ENV SIFDECODE
        ${CUTEst_DIR}/../SIFDecode
        ${CUTEst_DIR}/../sifdecode
)

find_program(SIFDecoder_EXE sifdecoder
    HINTS
        ${SIFDecoder_DIR}
    PATH_SUFFIXES
        bin
)

find_path(ARCHDefs_DIR compiler.${CUTEst_MYARCH}
    HINTS
        ${ARCHDEFS}
        ENV ARCHDEFS
        ${CUTEst_DIR}/../ARCHDefs
        ${CUTEst_DIR}/../archdefs
)

mark_as_advanced(SIFDecoder_DIR SIFDecoder_EXE ARCHDefs_DIR)
find_package_handle_standard_args(SIFDecoder 
    REQUIRED_VARS
        SIFDecoder_DIR
        SIFDecoder_EXE
        ARCHDefs_DIR
    NAME_MISMATCHED
)

# MASTSIF problems

find_path(MASTSIF_DIR mastsif.html
    HINTS
        ENV MASTSIF
        ${MASTSIF}
        ${CUTEst_DIR}/../MASTSIF
        ${CUTEst_DIR}/../mastsif
        ${CUTEst_DIR}/../SIF
        ${CUTEst_DIR}/../sif
)
mark_as_advanced(MASTSIF_DIR)
find_package_handle_standard_args(MASTSIF
    REQUIRED_VARS
        MASTSIF_DIR
    NAME_MISMATCHED
)
