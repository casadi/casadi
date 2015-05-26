message(STATUS "Looking for mosek")

find_path(MOSEK_INCLUDE_DIR
mosek.h
~/mosek/7/tools/platform/linux64x86/h
$ENV{MOSEK_INCLUDE_DIR}
)

find_library(MOSEK_LIBRARY
  NAMES mosek64
  PATHS ~/mosek/7/tools/platform/linux64x86/bin
  ${MOSEK_INCLUDE_DIR}/../bin ~/local/lib /usr/local/lib)

if(MOSEK_LIBRARY)
  set(FOUND_MOSEK_LIBS TRUE)
  set(MOSEK_LIBRARIES ${MOSEK_LIBRARY})
endif()

message("debug ${MOSEK_INCLUDE_DIR}")

set(MOSEK_FOUND TRUE)
