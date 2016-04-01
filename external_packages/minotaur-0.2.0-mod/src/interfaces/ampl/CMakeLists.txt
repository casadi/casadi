add_definitions(-DUSE_MINOTAUR_AMPL_INTERFACE)
include_directories("${PROJECT_BINARY_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/engines")
include_directories("${ASL_INC_DIR_F}")

set (ASL_LIB_SOURCES
  AMPLHessian.cpp
  AMPLInterface.cpp
  AMPLJacobian.cpp
  AMPLNonlinearFunction.cpp
)
set (ASL_LIB_HEADERS
  AMPLHessian.h
  AMPLInterface.h
  AMPLJacobian.h
  AMPLNonlinearFunction.cpp
)

add_library(mntrampl ${ASL_LIB_SOURCES})

# install the library at the user specified directory
if (BUILD_SHARED_LIBS)
  install(TARGETS mntrampl LIBRARY DESTINATION lib)
else()
  install(TARGETS mntrampl ARCHIVE DESTINATION lib)
endif()

# install the headers at the user specified directory
install(FILES ${ASL_LIB_HEADERS} DESTINATION include/minotaur)

# install the binary at the user specified directory
# install(TARGETS minotaurampl RUNTIME DESTINATION bin)

