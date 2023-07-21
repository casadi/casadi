set(CUTEST_MYARCH "pc64.lnx.gfo" CACHE STRING "CUTEst architecture" )

# CUTEst itself

find_path(CUTEST_INCLUDE_DIR cutest.h
    HINTS
        ${CUTEST}
        $ENV{CUTEST}
)
set(CUTEST_DIR "${CUTEST_INCLUDE_DIR}/.." CACHE PATH "")

find_library(CUTEST_LIBRARY
    NAMES
        cutest
    HINTS
        ${CUTEST_DIR}
    PATH_SUFFIXES
        objects/${CUTEST_MYARCH}/double
)

mark_as_advanced(CUTEST_INCLUDE_DIR CUTEST_DIR CUTEST_LIBRARY)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUTEst
    REQUIRED_VARS
        CUTEST_INCLUDE_DIR
        CUTEST_LIBRARY
)

# CUTEst library

if (CUTEst_FOUND AND NOT TARGET CUTEst::headers)
    add_library(CUTEst::headers INTERFACE IMPORTED)
    target_include_directories(CUTEst::headers
        INTERFACE ${CUTEST_INCLUDE_DIR})
endif()

if (CUTEst_FOUND AND NOT TARGET CUTEst::cutest)
    add_library(CUTEst::cutest STATIC IMPORTED)
    set_target_properties(CUTEst::cutest PROPERTIES
        IMPORTED_LOCATION ${CUTEST_LIBRARY})
    target_link_options(CUTEst::cutest INTERFACE
        LINKER:-whole-archive ${CUTEST_LIBRARY} LINKER:-no-whole-archive)
    target_link_libraries(CUTEst::cutest INTERFACE CUTEst::headers)
endif()

# if (CUTEst_FOUND AND NOT TARGET CUTEst::objects)
#     file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CUTEst)
#     set(CUTEst_OBJECTS
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cutest.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/pname.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/probname.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/varnames.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/newthread.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/problem.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/fortran_ops.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/interface.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ccutest.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/timings.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/usetup.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/udimen.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/udimse.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/udimsh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/unames.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/uvartype.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ufn.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ugr.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/uofg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/udh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ugrdh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ush.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ushp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ueh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ugreh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ugrsh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/uhprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ushprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ubandh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ureport.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/uterminate.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csetup.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cdimen.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cdimse.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cdimsh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cdimsj.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cdimchp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cnames.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cvartype.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cfn.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cgr.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cofg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cofsg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ccfg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/clfg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ccfsg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ccifg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ccifsg.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cdh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cdhc.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/ceh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cgrdh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cifn.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cigr.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cisgr.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cidh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cshc.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cshp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cish.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cjprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cstats.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csgr.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csgreh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csgrsh.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csjprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/chprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/chcprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cshprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cshcprod.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cchprods.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csjp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csgrp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/csgrshp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cchprodsp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/creport.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/connames.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cterminate.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/lqp.o
#         ${CMAKE_CURRENT_BINARY_DIR}/CUTEst/cconst.o)
#     add_custom_command(OUTPUT ${CUTEst_OBJECTS}
#         COMMAND ${CMAKE_AR} x ${CUTEST_LIBRARY}
#         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CUTEst
#         MAIN_DEPENDENCY ${CUTEST_LIBRARY})
#     add_custom_target(cutest-objects-extracted DEPENDS "${CUTEst_OBJECTS}")
#     add_library(CUTEst::objects OBJECT IMPORTED)
#     add_dependencies(CUTEst::objects cutest-objects-extracted)
#     set_target_properties(CUTEst::objects PROPERTIES
#         LINKER_LANGUAGE Fortran
#         IMPORTED_OBJECTS "${CUTEst_OBJECTS}")
#     target_link_libraries(CUTEst::objects INTERFACE CUTEst::headers)
# endif()
