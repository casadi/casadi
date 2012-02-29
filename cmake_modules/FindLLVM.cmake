find_program(LLVM_CONFIG NAMES llvm-config DOC "llvm-config executable")

if (LLVM_CONFIG)
  message(STATUS "LLVM llvm-config found at: ${LLVM_CONFIG}")
else (LLVM_CONFIG)
  message(FATAL_ERROR "Could NOT find llvm-config executable")
endif (LLVM_CONFIG)

execute_process(
  COMMAND ${LLVM_CONFIG} --includedir core
  OUTPUT_VARIABLE LLVM_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REPLACE " " ";" LLVM_INCLUDE_DIR ${LLVM_INCLUDE_DIR})

execute_process(
  COMMAND ${LLVM_CONFIG} --cppflags core
  OUTPUT_VARIABLE LLVM_CFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${LLVM_CONFIG} --libfiles core
  OUTPUT_VARIABLE LLVM_LIBRARIES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REPLACE " " ";" LLVM_LIBRARIES ${LLVM_LIBRARIES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LLVM DEFAULT_MSG LLVM_LIBRARIES LLVM_INCLUDE_DIR)
