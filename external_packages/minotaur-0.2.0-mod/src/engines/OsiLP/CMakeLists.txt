add_definitions(-DUSE_OSILP)
include_directories("${PROJECT_BINARY_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/base")


if (OSI_INC_DIR_F)
  include_directories("${OSI_INC_DIR_F}")
endif (OSI_INC_DIR_F)
if (CLP_INC_DIR_F)
  include_directories("${CLP_INC_DIR_F}")
endif (CLP_INC_DIR_F)
if (CPX_INC_DIR_F)
  include_directories("${CPX_INC_DIR}")
endif (CPX_INC_DIR_F)

add_definitions(-DCOIN_BIG_INDEX=0) ## avoid warnings
set (OSILP_SOURCES
  OsiLPEngine.cpp 
)
set (OSILP_HEADERS
  OsiLPEngine.h 
)

add_library(mntrosilp ${OSILP_SOURCES})

if (BUILD_SHARED_LIBS)
install(TARGETS mntrosilp LIBRARY DESTINATION lib)
else (BUILD_SHARED_LIBS)
install(TARGETS mntrosilp ARCHIVE DESTINATION lib)
endif (BUILD_SHARED_LIBS)

install(FILES ${OSILP_HEADERS} DESTINATION include/minotaur)

