cmake_minimum_required(VERSION 3.10.2)

# Three-stage build, as swig/rust:
#   SWIG_EXPORT  -- run SWIG only; emit casadiDOTNET_wrap.cxx + casadi.cs into gen/
#   SWIG_IMPORT  -- compile the pre-generated wrapper against an already-built core
#   neither      -- run SWIG inline, then build
# The native library (casadi_dotnet) needs only a C++ compiler; CasADi.dll
# (the C# assembly) is built when the dotnet CLI is found.

set(_iface    ${PROJECT_SOURCE_DIR}/swig/casadi.i)
set(_wrap_dir ${CMAKE_CURRENT_SOURCE_DIR}/gen)

if(SWIG_IMPORT)
  message(STATUS "SWIG_IMPORT: compiling pre-generated .NET wrapper; SWIG not invoked")
  set(_out_dir ${_wrap_dir})
else()
  execute_process(COMMAND ${SWIG_EXECUTABLE} -help OUTPUT_VARIABLE _swig_help ERROR_VARIABLE _swig_help_err)
  string(FIND "${_swig_help}${_swig_help_err}" "-dotnet" _has_dotnet)
  if(_has_dotnet EQUAL -1)
    message(FATAL_ERROR "SWIG at ${SWIG_EXECUTABLE} lacks -dotnet. Use the casadi-patched SWIG (ghcr.io/casadi/ci-swig).")
  endif()
  if(SWIG_EXPORT)
    set(_out_dir ${_wrap_dir})
  else()
    set(_out_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()
  file(MAKE_DIRECTORY ${_out_dir})
endif()

set(_swig_cpp ${_out_dir}/casadiDOTNET_wrap.cxx)
set(_swig_cs  ${_out_dir}/casadi.cs)

if(NOT SWIG_IMPORT)
  configure_file(../swig_config.h.cmake "${CMAKE_CURRENT_BINARY_DIR}/swig_config.h" ESCAPE_QUOTES)
  add_custom_command(
    OUTPUT  ${_swig_cpp} ${_swig_cs}
    COMMAND ${SWIG_EXECUTABLE} -dotnet -c++ -namespace CasADi -moduleclass ca -dllname casadi_dotnet
            ${CMAKE_SWIG_FLAGS} -w524,503 ${MAYBE_WERROR}
            -I${CMAKE_CURRENT_BINARY_DIR}
            -I${PROJECT_SOURCE_DIR}
            -I${PROJECT_SOURCE_DIR}/swig
            -outdir ${_out_dir}
            -o ${_swig_cpp} ${_iface}
    DEPENDS ${_iface} ${CASADI_SWIG_TARGET_FILES}
    COMMENT "swig -dotnet -> casadiDOTNET_wrap.cxx + casadi.cs"
    VERBATIM
  )
  add_custom_target(dotnet_source DEPENDS ${_swig_cpp} ${_swig_cs})
  if(SWIG_EXPORT)
    return()
  endif()
endif()

# ---- native library: plain C ABI, no .NET headers needed ----
add_library(casadi_dotnet SHARED ${_swig_cpp})
if(NOT SWIG_IMPORT)
  add_dependencies(casadi_dotnet dotnet_source)
endif()
if(WIN32)
  set_target_properties(casadi_dotnet PROPERTIES PREFIX "" CXX_STANDARD 11)
else()
  set_target_properties(casadi_dotnet PROPERTIES PREFIX "lib" CXX_STANDARD 11)
endif()
set_target_properties(casadi_dotnet PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                                               RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(casadi_dotnet PRIVATE ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_compile_definitions(casadi_dotnet PRIVATE WITH_DEPRECATED_FEATURES)
if(MSVC)
  target_compile_options(casadi_dotnet PRIVATE "/bigobj")
elseif(MINGW)
  target_compile_options(casadi_dotnet PRIVATE "-Wa,-mbig-obj")
endif()
target_link_libraries(casadi_dotnet casadi)
install(TARGETS casadi_dotnet DESTINATION ${LIB_PREFIX})
install(FILES ${_swig_cs} DESTINATION ${LIB_PREFIX}/dotnet)

# ---- CasADi.dll: generated casadi.cs, built by the dotnet CLI ----
find_program(DOTNET_EXECUTABLE dotnet)
if(DOTNET_EXECUTABLE)
  set(_asm_dir ${CMAKE_CURRENT_BINARY_DIR}/assembly)
  set(_asm ${_asm_dir}/CasADi.dll)
  add_custom_command(
    OUTPUT ${_asm}
    COMMAND ${DOTNET_EXECUTABLE} build ${CMAKE_CURRENT_SOURCE_DIR}/CasADi.csproj -c Release -nologo
            -p:CasadiGenerated=${_swig_cs} -o ${_asm_dir}
            -p:BaseIntermediateOutputPath=${CMAKE_CURRENT_BINARY_DIR}/obj/
    DEPENDS ${_swig_cs} ${CMAKE_CURRENT_SOURCE_DIR}/CasADi.csproj
    COMMENT "dotnet build -> CasADi.dll"
    VERBATIM
  )
  add_custom_target(casadi_dotnet_assembly ALL DEPENDS ${_asm})
  add_dependencies(casadi_dotnet_assembly casadi_dotnet)
  install(FILES ${_asm} DESTINATION ${LIB_PREFIX})
else()
  message(STATUS "No dotnet CLI found: CasADi.dll is not built (only casadi_dotnet)")
endif()
