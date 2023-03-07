# This file was partly taken from Sarcasm/irony-mode and modified

# Return variables
set(CLANG_LIBRARIES)
set(CLANG_DEFINITIONS)
set(CLANG_CXX_FLAGS)

# most recent versions come first: http://llvm.org/apt/
set(CLANG_KNOWN_LLVM_VERSIONS 4.0.1 4.0.0
  3.7.0 3.7
  3.6.2 3.6.1 3.6.0 3.6
  3.5.2 3.5.1 3.5.0 3.5
  3.4.2 3.4.1 3.4
  3.3
  3.2
  3.1)

# List of likely locations of llvm
set(llvm_hints)
foreach (version ${CLANG_KNOWN_LLVM_VERSIONS})
  string(REPLACE "." "" undotted_version "${version}")
  if(APPLE)
    # LLVM Homebrew
    set(llvm_hints ${llvm_hints} "/usr/local/Cellar/llvm/${version}/bin")
    # LLVM MacPorts
    set(llvm_hints ${llvm_hints} "/opt/local/libexec/llvm-${version}/bin")
  elseif(UNIX)
    # FreeBSD ports versions
    set(llvm_hints ${llvm_hints} "/usr/local/llvm${undotted_version}/bin")
  endif()
endforeach()

# Locate the LLVM config script
find_program(CLANG_LLVM_CONFIG llvm-config HINTS $ENV{CLANG}/bin ${llvm_hints})

# LLVM version
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --version
                OUTPUT_VARIABLE CLANG_LLVM_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Found LLVM ${CLANG_LLVM_VERSION}")

# LLVM installation prefix
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --prefix
                OUTPUT_VARIABLE CLANG_LLVM_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# LLVM include directory
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --includedir
                OUTPUT_VARIABLE CLANG_LLVM_INCLUDE_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# LLVM library directory
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --libdir
                OUTPUT_VARIABLE CLANG_LLVM_LIB_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# Clang shares include directory with llvm
set(CLANG_INCLUDE_DIR ${CLANG_LLVM_INCLUDE_DIR})

# All clang libraries
set(CLANG_DEP clangFrontend clangDriver clangCodeGen clangRewriteFrontend clangSerialization
              clangParse clangSema clangAnalysis clangEdit clangAST clangLex clangBasic)

# Get libraries
foreach(D ${CLANG_DEP})
  find_library(CLANG_DEP_${D} ${D} HINTS ${CLANG_LLVM_LIB_DIR})
  set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_DEP_${D}})
endforeach()

# Get the LLVM libraries
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --libfiles
                OUTPUT_VARIABLE CLANG_LLVM_LIBRARIES
                OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(CLANG_LLVM_LIBRARIES)
set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_LLVM_LIBRARIES})

# Get system libs
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --system-libs
                OUTPUT_VARIABLE CLANG_LLVM_SYSTEM_LIBS
                OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(CLANG_LLVM_SYSTEM_LIBS)
foreach(D ${CLANG_LLVM_SYSTEM_LIBS})
  if(${D} STREQUAL "-lm")
    # Standard math
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} m)
  elseif(${D} STREQUAL "-lz")
    # LLVM needs to be linked with zlib
    find_package(ZLIB QUIET REQUIRED)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${ZLIB_LIBRARIES})
  elseif(${D} STREQUAL "-lcurses")
    if(APPLE)
      # ??
      find_library(CLANG_CURSES ncurses)
    else()
      # ??
      find_library(CLANG_CURSES curses)
    endif()
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CLANG_CURSES})
  elseif(${D} STREQUAL "-lpthread")
    find_package(Threads QUIET REQUIRED)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  else()
    message(WARNING "Ignoring LLVM dependency '${D}'")
  endif()
endforeach()

# C++ compiler flags
message(STATUS "CLANG_LLVM_CONFIG: ${CLANG_LLVM_CONFIG}")
execute_process(COMMAND ${CLANG_LLVM_CONFIG} --cxxflags
                OUTPUT_VARIABLE CLANG_LLVM_CXXFLAGS
                OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(CLANG_LLVM_CXXFLAGS)
foreach(D ${CLANG_LLVM_CXXFLAGS})
  if(${D} MATCHES "-D")
    set(CLANG_DEFINITIONS ${CLANG_DEFINITIONS} ${D})
  else()
    message(STATUS "Ignoring LLVM C++ flag '${D}'")
  endif()
endforeach()

if(MINGW)
  message("/usr/${PREFIX}/lib")
  find_library(IMAGEHLP imagehlp /usr/${PREFIX}/lib)
  if (IMAGEHLP)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${IMAGEHLP})
  endif()
  find_library(SHLWAPI libshlwapi /usr/${PREFIX}/lib)
  if (SHLWAPI)
    set(CLANG_LIBRARIES ${CLANG_LIBRARIES} ${SHLWAPI})
  endif()
  message("${CLANG_LIBRARIES}")
endif()


set(CLANG_FOUND TRUE)


add_library(clang::clang INTERFACE IMPORTED)
set_target_properties(clang::clang PROPERTIES
  INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${CLANG_INCLUDE_DIR}
)
target_link_libraries(clang::clang INTERFACE ${CLANG_LIBRARIES})

