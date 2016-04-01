add_definitions(-DUSE_IPOPT)
include_directories("${PROJECT_BINARY_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/base")
include_directories("${IPOPT_INC_DIR_F}")

set (IPOPT_SOURCES
  IpoptDummy.f
  IpoptEngine.cpp 
)
set (IPOPT_HEADERS
  IpoptEngine.h 
  IpoptEngineTnlp.h
)

add_library(mntripopt ${IPOPT_SOURCES})

if (BUILD_SHARED_LIBS)
install(TARGETS mntripopt LIBRARY DESTINATION lib)
else (BUILD_SHARED_LIBS)
install(TARGETS mntripopt ARCHIVE DESTINATION lib)
endif (BUILD_SHARED_LIBS)

install(FILES ${IPOPT_HEADERS} DESTINATION include/minotaur)

