##############################################################################
## Do not change this section. It should stay same for all binaries
##############################################################################
add_definitions(-DUSE_MINOTAUR_AMPL_INTERFACE)
include_directories("${PROJECT_BINARY_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/engines")
include_directories("${PROJECT_SOURCE_DIR}/src/interfaces/ampl")
include_directories("${ASL_INC_DIR_F}")

set (QPOASES_LIBS)
set (GFORTRAN_COMPILER)

## define where to search for external libraries. This path must be defined
## before add_executable is called.
if (LINK_BQPD)
  add_definitions(-DUSE_BQPD)
  include_directories("${PROJECT_SOURCE_DIR}/src/engines/Bqpd")
  link_directories(${BQPD_LIB_DIR_F})
endif()

if (LINK_FILTER_SQP)
  add_definitions(-DUSE_FilterSQP)
  include_directories("${PROJECT_SOURCE_DIR}/src/engines/FilterSQP")
  link_directories(${FILTER_LIB_DIR_F})
endif()

if (LINK_IPOPT)
  add_definitions(-DUSE_IPOPT)
  include_directories("${PROJECT_SOURCE_DIR}/src/engines/Ipopt")
  include_directories("${IPOPT_INC_DIR_F}")
  link_directories(${IPOPT_LIB_DIR_F})
endif()

if (LINK_OSI)
  add_definitions(-DUSE_OSILP)
  include_directories("${PROJECT_SOURCE_DIR}/src/engines/OsiLP")
  if (OSI_INC_DIR_F)
    include_directories("${OSI_INC_DIR_F}")
  endif()
  if (OSI_LIB_DIR_F)
    link_directories(${OSI_LIB_DIR_F})
  endif()
  if (LINK_CLP)
    link_directories(${CLP_LIB_DIR_F})
  endif()
  if (LINK_CPX)
    link_directories(${CPX_LIB_DIR})
  endif()
endif()

if (LINK_QPOASES)
  add_definitions(-DUSE_QPOASES)
  include_directories("${PROJECT_SOURCE_DIR}/src/engines/qpOASES")
  include_directories("${QPOASES_INC_DIR}")
  set (QPOASES_LIBS mntrqpoases)
endif()
  
if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  if (LINK_IPOPT OR LINK_FILTER_SQP OR LINK_BQPD)
    set (GFORTRAN_COMPILER gfortran)
  endif()
endif()

set (ALL_EXEC_LIBS 
	"mntrampl" 
	"mntrengfac" 
	"${IPOPT_LIBS}"
        "${OSI_LIBS}"
        "${FILTER_SQP_LIBS}"
	"${BQPD_LIBS}"
        "${QPOASES_LIBS}"
        "minotaur"
	"${ASL_LIBS}"
        "lapack"
	"blas"
        "${GFORTRAN_COMPILER}"
        "${MNTR_EXTRA_LIBS}"
	)
message(STATUS ${MSG_HEAD} "ALL_EXEC_LIBS =  ${ALL_EXEC_LIBS}")
list(REVERSE ALL_EXEC_LIBS)
list(REMOVE_DUPLICATES ALL_EXEC_LIBS)
list(REVERSE ALL_EXEC_LIBS)
message(STATUS ${MSG_HEAD} "ALL_EXEC_LIBS after removing duplicates =  ${ALL_EXEC_LIBS}")


##############################################################################
## Add lines specific to your binaries in this section. 
## Use the lines meant for bnb as a template.
##############################################################################
set (BNB_SOURCES
  Bnb.cpp 
)

add_executable(bnb ${BNB_SOURCES})
target_link_libraries(bnb ${ALL_EXEC_LIBS})

# This will install the binary in bin directory.
install(TARGETS bnb RUNTIME DESTINATION bin)

##############################################################################
## Add lines specific to your binaries in this section.
## Use the lines meant for bnb as a template.
##############################################################################
#set (MSBNB_SOURCES
#  MsBnb.cpp
#)
#
#add_executable(msbnb ${MSBNB_SOURCES})
#target_link_libraries(msbnb ${ALL_EXEC_LIBS})
#
## This will install the binary in bin directory.
#install(TARGETS msbnb RUNTIME DESTINATION bin)


##############################################################################
## Add lines specific to your binaries in this section. 
## Use the lines meant for bnb as a template.
##############################################################################
set (QPD_SOURCES
  QPDive.cpp 
)

add_executable(qpd ${QPD_SOURCES})
target_link_libraries(qpd ${ALL_EXEC_LIBS})

# This will install the binary in bin directory.
install(TARGETS qpd RUNTIME DESTINATION bin)

##############################################################################
## Add lines specific to your binaries in this section. 
## Use the lines meant for bnb as a template.
##############################################################################
set (QG_SOURCES
  qg.cpp 
)

add_executable(qg ${QG_SOURCES})
target_link_libraries(qg ${ALL_EXEC_LIBS})

# This will install the binary in bin directory.
install(TARGETS qg RUNTIME DESTINATION bin)



##############################################################################
## Add lines specific to your binaries in this section. 
## Use the lines meant for bnb as a template.
##############################################################################
#set (POLYNOMIAL_SOURCES
#  polynomial-solver.cpp 
#)
#
#add_executable(polly ${POLYNOMIAL_SOURCES})
#target_link_libraries(polly ${ALL_EXEC_LIBS})
#
## This will install the binary in bin directory.
#install(TARGETS polly RUNTIME DESTINATION bin)


##############################################################################
## Add lines specific to your binaries in this section. 
## Use the lines meant for bnb as a template.
##############################################################################
set (GLOB_SOURCES
  Glob.cpp 
)

add_executable(glob ${GLOB_SOURCES})
target_link_libraries(glob ${ALL_EXEC_LIBS})

# This will install the binary in bin directory.
install(TARGETS glob RUNTIME DESTINATION bin)

# Serdar added for branch and cut method. 
##############################################################################
## Add lines specific to your binaries in this section. 
## Use the lines meant for bnb as a template.
##############################################################################
#set (BNC_SOURCES
#  Bnc.cpp 
#)
#
#add_executable(bnc ${BNC_SOURCES})
#target_link_libraries(bnc ${ALL_EXEC_LIBS})
#
## This will install the binary in bin directory.
#install(TARGETS bnc RUNTIME DESTINATION bin)
