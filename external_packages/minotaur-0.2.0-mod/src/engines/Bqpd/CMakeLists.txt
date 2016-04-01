add_definitions(-DUSE_BQPD)
include_directories("${PROJECT_BINARY_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/base")

set (BQPD_SOURCES
  BqpdEngine.cpp 
  BqpdAux.f 
  BqpdAux2.f 
)
set (BQPD_HEADERS
  BqpdEngine.h 
)


add_library(mntrbqpd ${BQPD_SOURCES})

if (BUILD_SHARED_LIBS)
install(TARGETS mntrbqpd LIBRARY DESTINATION lib)
else (BUILD_SHARED_LIBS)
install(TARGETS mntrbqpd ARCHIVE DESTINATION lib)
endif (BUILD_SHARED_LIBS)

install(FILES ${BQPD_HEADERS} DESTINATION include/minotaur)

