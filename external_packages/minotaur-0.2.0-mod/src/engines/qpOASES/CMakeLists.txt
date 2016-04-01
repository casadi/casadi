add_definitions(-DUSE_QPOASES)
include_directories("${PROJECT_BINARY_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/base")

include_directories("${QPOASES_INC_DIR}")

set (QPOASES_SOURCES
  qpOASESEngine.cpp
)
set (QPOASES_HEADERS
  qpOASESEngine.h
)


add_library(mntrqpoases ${QPOASES_SOURCES})

target_link_libraries (mntrqpoases ${QPOASES_LIB_DIR}/libqpOASES.a)

if (BUILD_SHARED_LIBS)
install(TARGETS mntrqpoases LIBRARY DESTINATION lib)
else (BUILD_SHARED_LIBS)
install(TARGETS mntrqpoases ARCHIVE DESTINATION lib)
endif (BUILD_SHARED_LIBS)

install(FILES ${QPOASES_HEADERS} DESTINATION include/minotaur)

