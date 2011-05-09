find_package(LAPACK)
find_package(Threads)

FIND_PATH(IPOPT_INCLUDE_DIR 
coin/IpIpoptApplication.hpp
HINTS $ENV{IPOPT}/include
)

IF (IPOPT_INCLUDE_DIR)
   SET(IPOPT_FOUND_INCLUDE TRUE)
   MESSAGE(STATUS "Found IPOPT include dir: ${IPOPT_INCLUDE_DIR}")
ELSE (IPOPT_INCLUDE_DIR)
   MESSAGE(STATUS "Could not find IPOPT include dir")
ENDIF (IPOPT_INCLUDE_DIR)

FIND_LIBRARY(IPOPT_LIBRARY 
ipopt HINTS $ENV{IPOPT}/coin/lib/ $ENV{IPOPT}/coin/lib/ThirdParty/ $ENV{IPOPT}/lib/ $ENV{IPOPT}/lib/coin/  /usr/local/lib/coin /usr/local/lib/coin/ThirdParty )

FIND_LIBRARY(HSL_LIBRARY 
coinhsl HINTS  $ENV{IPOPT}/coin/lib/ $ENV{IPOPT}/coin/lib/ThirdParty/ $ENV{IPOPT}/lib/ $ENV{IPOPT}/lib/coin/  /usr/local/lib/coin /usr/local/lib/coin/ThirdParty)


FIND_LIBRARY(MUMPS_LIBRARY_492 
dmumps_seq-4.9.2)

MESSAGE(STATUS "${PTHREAD_LIBRARIES}")

IF (IPOPT_LIBRARY)
   SET(IPOPT_FOUND_LIBS TRUE)
   SET(IPOPT_LIBRARIES ${IPOPT_LIBRARY} ${IPOPT_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${LAPACK_LIBRARIES} ${EXTRA_LIBRARIES} ${CMAKE_DL_LIBS})

  IF(MUMPS_LIBRARY_492)
    SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${MUMPS_LIBRARY_492})
  ELSE(MUMPS_LIBRARY_492)
    MESSAGE(STATUS "MUMPS warning: If you have installed IPOPT through a package manager and receive an error ImportError: /usr/lib/libipopt.so.0: undefined symbol: MPI_Init when using IPOPT, you should install libmumps-seq-4.9.2")
  ENDIF(MUMPS_LIBRARY_492)
  
  IF(HSL_LIBRARY)
    SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${HSL_LIBRARY})
  ENDIF(HSL_LIBRARY)

   MESSAGE(STATUS "Found Ipopt libs: ${IPOPT_LIBRARIES}")
ELSE (IPOPT_LIBRARY)
   MESSAGE(STATUS "Could not find Ipopt libs")
ENDIF (IPOPT_LIBRARY)

IF(IPOPT_FOUND_INCLUDE AND IPOPT_FOUND_LIBS)
  SET(IPOPT_FOUND TRUE)
ENDIF(IPOPT_FOUND_INCLUDE AND IPOPT_FOUND_LIBS)





