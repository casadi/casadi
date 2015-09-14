if (OLD_LLVM)
# Legacy version, to be removed
message(STATUS "Looking for clang 3.4.2")

set(LLVM_INSTALL_PREFIX $ENV{CLANG})

set(LLVM_DEP LLVMJIT LLVMInterpreter LLVMX86CodeGen LLVMBitWriter LLVMIRReader LLVMipo LLVMLinker LLVMRuntimeDyld LLVMExecutionEngine LLVMAsmPrinter LLVMSelectionDAG LLVMX86Desc LLVMAsmParser LLVMBitReader LLVMVectorize LLVMMCParser LLVMCodeGen LLVMX86AsmPrinter LLVMX86Info LLVMObjCARCOpts LLVMScalarOpts LLVMX86Utils LLVMInstCombine LLVMTransformUtils LLVMipa LLVMAnalysis LLVMTarget LLVMCore LLVMMC LLVMObject LLVMSupport LLVMBitWriter LLVMIRReader LLVMipo LLVMLinker LLVMVectorize LLVMInstrumentation LLVMBitReader LLVMOption LLVMX86CodeGen LLVMAsmParser LLVMAArch64AsmParser LLVMAArch64Disassembler LLVMARMCodeGen LLVMARMAsmParser LLVMARMDisassembler LLVMCppBackendCodeGen LLVMHexagonCodeGen LLVMMipsCodeGen LLVMMipsAsmParser LLVMMipsDisassembler LLVMMSP430CodeGen LLVMNVPTXCodeGen LLVMPowerPCCodeGen LLVMPowerPCAsmParser LLVMR600CodeGen LLVMSparcCodeGen LLVMSystemZCodeGen LLVMSystemZAsmParser LLVMSystemZDisassembler LLVMX86AsmParser LLVMX86Disassembler LLVMX86Desc LLVMX86AsmPrinter LLVMX86Utils LLVMX86Info LLVMXCoreCodeGen LLVMXCoreDisassembler LLVMAArch64CodeGen LLVMAsmPrinter LLVMMCParser LLVMSelectionDAG LLVMCodeGen LLVMObjCARCOpts LLVMScalarOpts LLVMInstCombine LLVMTransformUtils LLVMipa LLVMAnalysis LLVMARMDesc LLVMCppBackendInfo LLVMHexagonAsmPrinter LLVMMipsDesc LLVMMSP430Desc LLVMNVPTXDesc LLVMPowerPCDesc LLVMR600Desc LLVMSparcDesc LLVMSystemZDesc LLVMXCoreDesc LLVMAArch64Desc LLVMARMAsmPrinter LLVMARMInfo LLVMHexagonDesc LLVMMipsAsmPrinter LLVMMipsInfo LLVMMSP430AsmPrinter LLVMMSP430Info LLVMNVPTXAsmPrinter LLVMNVPTXInfo LLVMPowerPCAsmPrinter LLVMPowerPCInfo LLVMR600AsmPrinter LLVMR600Info LLVMSparcInfo LLVMSystemZAsmPrinter LLVMSystemZInfo LLVMXCoreAsmPrinter LLVMXCoreInfo LLVMAArch64AsmPrinter LLVMAArch64Info LLVMTarget LLVMHexagonInfo LLVMAArch64Utils LLVMCore LLVMMC LLVMObject LLVMSupport)

set(LLVM_DEFINITIONS "-DCLANG_ENABLE_OBJC_REWRITER" "-DGTEST_HAS_RTTI=0" "-D_GNU_SOURCE" "-D__STDC_CONSTANT_MACROS" "-D__STDC_FORMAT_MACROS" "-D__STDC_LIMIT_MACROS")
set(CLANG_CXX_FLAGS "-fPIC -fvisibility-inlines-hidden -ffunction-sections -fdata-sections -fno-common -fno-strict-aliasing")

else (OLD_LLVM)

# LLVM needs to be linked with zlib (?)
find_package(ZLIB REQUIRED)
message(STATUS "ZLIB libraries ${ZLIB_LIBRARIES}")

message(STATUS "Looking clang 3.5 or newer")

find_package(LLVM REQUIRED CONFIG)

message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")
message(STATUS "LLVM install prefix ${LLVM_INSTALL_PREFIX}")
message(STATUS "LLVM libraries ${LLVM_INSTALL_PREFIX}")

set(LLVM_DEP ${LLVM_AVAILABLE_LIBS})
set(CLANG_CXX_FLAGS)

endif (OLD_LLVM)

# Clang shares include directory with llvm
set(CLANG_INCLUDE_DIR ${LLVM_INSTALL_PREFIX}/include)

# All clang libraries
set(CLANG_DEP clangFrontend clangDriver clangCodeGen clangRewriteFrontend clangFrontend clangSerialization clangParse clangSema clangAnalysis clangEdit clangAST clangLex clangBasic ${LLVM_DEP})

# Get libraries
set(CLANG_LIBRARIES)
foreach(D ${CLANG_DEP})
  find_library(CLANG_DEP_${D} ${D} ${LLVM_INSTALL_PREFIX}/lib)
  set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_DEP_${D}})
endforeach()

if (ZLIB_LIBRARIES)
  set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${ZLIB_LIBRARIES})
endif()

if(MINGW)
  message("/usr/${PREFIX}/lib")
  find_library(CLANG_EXTRA_LIBS imagehlp /usr/${PREFIX}/lib)
  set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_EXTRA_LIBS})
  message("${CLANG_LIBRARIES}")
endif()

if(APPLE)
  find_library(CLANG_EXTRA_LIBS ncurses)
  set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_EXTRA_LIBS})
endif()

set(CLANG_FOUND TRUE)

