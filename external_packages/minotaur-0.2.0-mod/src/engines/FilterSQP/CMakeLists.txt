add_definitions(-DUSE_FILTERSQP)
include_directories("${PROJECT_BINARY_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/base")

set (FILTERSQP_SOURCES
  FilterSQPEngine.cpp 
)
set (FILTERSQP_HEADERS
  FilterSQPEngine.h 
)

add_library(mntrfiltersqp ${FILTERSQP_SOURCES})

if (BUILD_SHARED_LIBS)
install(TARGETS mntrfiltersqp LIBRARY DESTINATION lib)
else (BUILD_SHARED_LIBS)
install(TARGETS mntrfiltersqp ARCHIVE DESTINATION lib)
endif (BUILD_SHARED_LIBS)


install(FILES ${FILTERSQP_HEADERS} DESTINATION include/minotaur)

