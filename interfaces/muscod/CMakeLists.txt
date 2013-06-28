cmake_minimum_required (VERSION 2.6)

set(ROOT /home/janderss/src/MUSCOD)
set(LACPATH   ${ROOT}/MUSCOD-II/LIBLAC)
set(MC2PATH   ${ROOT}/MUSCOD-II/MC2)
set(PGPLOTPATH ${ROOT}/MUSCOD-II/PGPLOT)
set(PROBLEMPATH   ${ROOT}/MUSCOD-II/MC2_TEST)

include_directories(${MC2PATH}/INCLUDE ${LACPATH}/INCLUDE ${PROBLEMPATH}/SRC)

link_directories(${PROBLEMPATH} ${PROBLEMPATH}/SRC ${MC2PATH}/IND ${MC2PATH}/QPS ${MC2PATH}/MINLP ${LACPATH} ${PGPLOTPATH} /usr/X11R6/lib)

set(CCUFLAGS "-g -DNDEBUG -DCSTATS -DMSPLOT -DLSQ -DCLSQ -DPRSQP")
set(FCUFLAGS "-g")
set(CXXCUFLAGS ${CCUFLAGS})

# Defines call of make in subdirectories. Use option -s for silent mode
set(CCSFLAGS "-DLINUX -Wall -D__x86__ -D__linux__ -D__OSVERSION__=2")
set(CXXCSFLAGS "-DLINUX -DMALLOC -Wno-deprecated -Wno-unused ${CCSFLAGS}")
set(CPFLAGS "-P -C -DLINUX -traditional-cpp")
set(FCSFLAGS "-Wall -fno-second-underscore -fno-automatic -finit-local-zero")

set(CFLAGS    "${CCSFLAGS} ${CCUFLAGS} ${IFLAGS}")
set(FFLAGS    "${FCSFLAGS} ${FCUFLAGS}")
set(CXXFLAGS   "${CXXCSFLAGS} ${CCUFLAGS} ${IFLAGS}")

set(CMAKE_C_FLAGS ${CFLAGS})
set(CMAKE_CXX_FLAGS ${CXXFLAGS})


set(MUSCOD_INTERFACE_SRCS
  muscod_interface.hpp
  muscod_interface.cpp
  muscod_internal.hpp
  muscod_internal.cpp
  muscod_function.hpp
  muscod_function.cpp
  muscod_aux.hpp
  muscod_aux.cpp
)

if(ENABLE_STATIC)
add_library(casadi_muscod_interface STATIC ${MUSCOD_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_muscod_interface SHARED ${MUSCOD_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_muscod_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./
   DESTINATION include/casadi/interfaces/muscod
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

