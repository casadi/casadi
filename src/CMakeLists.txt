find_package(Eigen3 REQUIRED)

function(alpaqa_configure_visibility target)
    set_target_properties(${target} PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                               C_VISIBILITY_PRESET "hidden"
                                               VISIBILITY_INLINES_HIDDEN true)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        target_link_options(${target} PRIVATE
            $<$<LINK_LANGUAGE:C,CXX>:LINKER:--exclude-libs,ALL>)
    endif()
endfunction()

option(ALPAQA_WITH_OCP 
    "Enable solvers tailored for optimal control problems" On)

add_library(alpaqa
    "alpaqa/src/util/type-erasure.cpp"
    "alpaqa/src/util/demangled-typename.cpp"
    "alpaqa/src/util/print.cpp"
    "alpaqa/src/util/io/csv.cpp"
    "alpaqa/src/accelerators/lbfgs.cpp"
    "alpaqa/src/problem/problem-counters.cpp"
    "alpaqa/src/problem/ocproblem-counters.cpp"
    "alpaqa/src/problem/type-erased-problem.cpp"
    "alpaqa/src/outer/alm.cpp"
    "alpaqa/src/outer/internal/alm-helpers.cpp"
    "alpaqa/src/inner/panoc.cpp"
    "alpaqa/src/inner/zerofpr.cpp"
    "alpaqa/src/inner/pantr.cpp"
    "alpaqa/src/panoc-alm.cpp"
    "alpaqa/src/panoc-anderson-alm.cpp"
    "alpaqa/src/structured-panoc-alm.cpp"
    "alpaqa/src/structured-zerofpr-alm.cpp"
    "alpaqa/src/zerofpr-alm.cpp"
    "alpaqa/src/zerofpr-anderson-alm.cpp"
    "alpaqa/src/newton-tr-pantr-alm.cpp"
    "alpaqa/src/inner/internal/solverstatus.cpp"
    "alpaqa/src/inner/directions/panoc/structured-lbfgs.cpp"
    "alpaqa/src/inner/directions/panoc/structured-newton.cpp"
    "alpaqa/src/inner/internal/panoc-helpers.cpp"
)
if (ALPAQA_WITH_OCP)
    target_sources(alpaqa PRIVATE
        "alpaqa/src/inner/panoc-ocp.cpp"
        "alpaqa/src/problem/ocproblem.cpp"
    )
endif()
if (ALPAQA_WITH_DRIVERS)
    target_sources(alpaqa PRIVATE "alpaqa/src/params/params.cpp")
endif()
target_compile_features(alpaqa PUBLIC cxx_std_20)
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/alpaqa/include>
    $<INSTALL_INTERFACE:include>)
target_compile_definitions(alpaqa PUBLIC
    $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
    $<$<CONFIG:Debug>:EIGEN_RUNTIME_NO_MALLOC>)
target_link_libraries(alpaqa PUBLIC Eigen3::Eigen)
target_link_libraries(alpaqa PRIVATE warnings)
alpaqa_configure_visibility(alpaqa)
target_compile_definitions(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:ALPAQA_WITH_QUAD_PRECISION>)
target_compile_definitions(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_OCP}>:ALPAQA_WITH_OCP>)
target_link_libraries(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:quadmath>)
set_target_properties(alpaqa PROPERTIES SOVERSION ${PROJECT_VERSION})
add_library(alpaqa::alpaqa ALIAS alpaqa)
list(APPEND ALPAQA_INSTALL_TARGETS alpaqa)

# DLL import/export
include(GenerateExportHeader)
generate_export_header(alpaqa
    EXPORT_FILE_NAME export/alpaqa/export.h)
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)

# CasADi
if (ALPAQA_WITH_CASADI)
    find_package(casadi REQUIRED)
    find_package(Threads REQUIRED)
    # Normal NLPs
    add_library(casadi-loader 
        "interop/casadi/src/CasADiProblem.cpp"
        "interop/casadi/include/alpaqa/casadi/CasADiFunctionWrapper.hpp")
    target_include_directories(alpaqa PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/casadi/include>
        $<INSTALL_INTERFACE:include>)
    set_property(TARGET casadi-loader
        PROPERTY OUTPUT_NAME alpaqa-casadi-loader)
    target_link_libraries(casadi-loader
        PUBLIC  alpaqa::alpaqa
        PRIVATE casadi Threads::Threads alpaqa::warnings)
    alpaqa_configure_visibility(casadi-loader)
    set_target_properties(casadi-loader PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::casadi-loader ALIAS casadi-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS casadi-loader)
    generate_export_header(casadi-loader
        EXPORT_FILE_NAME export/alpaqa/casadi-loader-export.h)
    target_include_directories(casadi-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)

    message(STATUS "Compiling WITH CasADi support")
else()
    message(STATUS "Compiling WITHOUT CasADi support")
endif()

set(ALPAQA_WITH_CASADI_OCP_SUPPORTED ${ALPAQA_WITH_OCP})
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16.0.0)
    set(ALPAQA_WITH_CASADI_OCP_SUPPORTED Off)
endif()

include(CMakeDependentOption)
cmake_dependent_option(ALPAQA_WITH_CASADI_OCP
    "Build the CasADi loader for optimal control problems"
    ${ALPAQA_WITH_CASADI_OCP_SUPPORTED} "ALPAQA_WITH_CASADI" OFF)

# Optimal control problems
if (ALPAQA_WITH_OCP AND ALPAQA_WITH_CASADI_OCP)
    add_library(casadi-ocp-loader "interop/casadi/src/CasADiControlProblem.cpp")
    target_include_directories(casadi-ocp-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/casadi/include>
        $<INSTALL_INTERFACE:include>)
    set_property(TARGET casadi-ocp-loader
        PROPERTY OUTPUT_NAME alpaqa-casadi-ocp-loader)
    target_link_libraries(casadi-ocp-loader
        PUBLIC  alpaqa::alpaqa
        PRIVATE casadi Threads::Threads alpaqa::warnings)
    alpaqa_configure_visibility(casadi-ocp-loader)
    set_target_properties(casadi-ocp-loader PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::casadi-ocp-loader ALIAS casadi-ocp-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS casadi-ocp-loader)
    generate_export_header(casadi-ocp-loader
        EXPORT_FILE_NAME export/alpaqa/casadi-ocp-loader-export.h)
    target_include_directories(casadi-ocp-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

if (CMAKE_DL_LIBS)
    add_library(dl-api INTERFACE)
    target_include_directories(dl-api INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/dl-api/include>
        $<INSTALL_INTERFACE:include>)
    add_library(alpaqa::dl-api ALIAS dl-api)
    list(APPEND ALPAQA_INSTALL_TARGETS dl-api)
    add_library(dl-loader "interop/dl/src/dl-problem.cpp")
    target_include_directories(dl-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/dl/include>
        $<INSTALL_INTERFACE:include>)
    set_property(TARGET dl-loader
        PROPERTY OUTPUT_NAME alpaqa-dl-loader)
    target_link_libraries(dl-loader
        PUBLIC alpaqa::alpaqa alpaqa::dl-api
        PRIVATE ${CMAKE_DL_LIBS} alpaqa::warnings)
    set_target_properties(dl-loader PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::dl-loader ALIAS dl-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS dl-loader)
    message(STATUS "Compiling WITH dl support")
else()
    message(STATUS "Compiling WITHOUT dl support")
endif()

if (ALPAQA_WITH_CUTEST)
    find_package(CUTEst REQUIRED)
    add_library(cutest-interface "interop/cutest/src/cutest-loader.cpp")
    target_include_directories(cutest-interface PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/cutest/include>
        $<INSTALL_INTERFACE:include>)
    target_link_libraries(cutest-interface
        PUBLIC alpaqa::alpaqa
        PRIVATE CUTEst::headers ${CMAKE_DL_LIBS} alpaqa::warnings)
    set_target_properties(cutest-interface PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::cutest-interface ALIAS cutest-interface)
    list(APPEND ALPAQA_INSTALL_TARGETS cutest-interface)
    message(STATUS "Compiling WITH CUTEst support")
else()
    message(STATUS "Compiling WITHOUT CUTEst support")
endif()

find_package(Ipopt)

if (TARGET COIN::Ipopt)
    message(STATUS "Compiling Ipopt problem adapter")
    add_library(ipopt-adapter
        "interop/ipopt/src/ipopt-adapter.cpp"
        "interop/ipopt/src/ipopt-params.cpp"
    )
    target_link_libraries(ipopt-adapter 
        PUBLIC COIN::Ipopt alpaqa::alpaqa PRIVATE alpaqa::warnings)
    target_link_directories(ipopt-adapter PUBLIC
        $<BUILD_INTERFACE:${IPOPT_LIBRARY_DIRS}>)
    target_include_directories(ipopt-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/ipopt/include>
        $<INSTALL_INTERFACE:include>)
    set_property(TARGET ipopt-adapter
        PROPERTY OUTPUT_NAME alpaqa-ipopt-adapter)
    alpaqa_configure_visibility(ipopt-adapter)
    set_target_properties(ipopt-adapter PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::ipopt-adapter ALIAS ipopt-adapter)
    list(APPEND ALPAQA_INSTALL_TARGETS ipopt-adapter)
    generate_export_header(ipopt-adapter
        EXPORT_FILE_NAME export/alpaqa/ipopt-adapter-export.h)
    target_include_directories(ipopt-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

find_package(lbfgspp QUIET)
if (TARGET lbfgspp)
    message(STATUS "Compiling LBFGS++ solver adapter")
    add_library(lbfgspp-adapter "interop/lbfgspp/src/lbfgsb-adapter.cpp")
    target_link_libraries(lbfgspp-adapter 
        PUBLIC lbfgspp alpaqa::alpaqa PRIVATE alpaqa::warnings)
    target_include_directories(lbfgspp-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/lbfgspp/include>
        $<INSTALL_INTERFACE:include>)
    set_property(TARGET lbfgspp-adapter
        PROPERTY OUTPUT_NAME alpaqa-lbfgspp-adapter)
    alpaqa_configure_visibility(lbfgspp-adapter)
    set_target_properties(lbfgspp-adapter PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::lbfgspp-adapter ALIAS lbfgspp-adapter)
    list(APPEND ALPAQA_INSTALL_TARGETS lbfgspp-adapter)
    generate_export_header(lbfgspp-adapter
        EXPORT_FILE_NAME export/alpaqa/lbfgspp-adapter-export.h)
    target_include_directories(lbfgspp-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

cmake_dependent_option(ALPAQA_WITH_LBFGSB
    "Include the L-BFGS-B solver" ON ALPAQA_HAVE_FORTRAN OFF)

if (ALPAQA_WITH_LBFGSB)
    message(STATUS "Compiling L-BFGS-B solver adapter")
    add_library(lbfgsb-fortran STATIC
        "thirdparty/lbfgsb/Lbfgsb.3.0/lbfgsb.f"
        "thirdparty/lbfgsb/Lbfgsb.3.0/timer.f"
        "thirdparty/lbfgsb/Lbfgsb.3.0/blas.f"
        "thirdparty/lbfgsb/Lbfgsb.3.0/linpack.f"
    )
    set_target_properties(lbfgsb-fortran PROPERTIES LINKER_LANGUAGE Fortran)
    list(APPEND ALPAQA_INSTALL_TARGETS lbfgsb-fortran)
    add_library(lbfgsb-adapter
        "interop/lbfgsb/src/lbfgsb-adapter.cpp"
        "interop/lbfgsb/src/lbfgsb-params.cpp"
        "interop/lbfgsb/src/lbfgsb_c.f90"
    )
    set_property(SOURCE "interop/lbfgsb/src/lbfgsb_c.f90"
                 PROPERTY Fortran_PREPROCESS On)
    target_link_libraries(lbfgsb-adapter 
        PUBLIC lbfgsb-fortran alpaqa::alpaqa PRIVATE alpaqa::warnings)
    target_include_directories(lbfgsb-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/lbfgsb/include>
        $<INSTALL_INTERFACE:include>)
    set_property(TARGET lbfgsb-adapter
        PROPERTY OUTPUT_NAME alpaqa-lbfgsb-adapter)
    alpaqa_configure_visibility(lbfgsb-adapter)
    set_target_properties(lbfgsb-adapter PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::lbfgsb-adapter ALIAS lbfgsb-adapter)
    list(APPEND ALPAQA_INSTALL_TARGETS lbfgsb-adapter)
    generate_export_header(lbfgsb-adapter
        EXPORT_FILE_NAME export/alpaqa/lbfgsb-adapter-export.h)
    target_include_directories(lbfgsb-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

# Drivers
if (ALPAQA_WITH_DRIVERS)
    add_executable(driver
        "alpaqa/src/driver/alpaqa-driver.cpp"
        "alpaqa/src/driver/panoc-driver.cpp"
        "alpaqa/src/driver/pantr-driver.cpp"
        "alpaqa/src/driver/ipopt-driver.cpp"
        "alpaqa/src/driver/lbfgsb-driver.cpp"
        "alpaqa/src/driver/problem.cpp"
    )
    set_target_properties(driver PROPERTIES
        OUTPUT_NAME "alpaqa-driver"
        RELEASE_POSTFIX ""
        DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}"
        RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX ""
    )
    target_link_libraries(driver
        PRIVATE alpaqa::alpaqa alpaqa::warnings alpaqa::dl-loader)
    add_executable(alpaqa::driver ALIAS driver)
    list(APPEND ALPAQA_INSTALL_EXE driver)
    if (TARGET alpaqa::casadi-loader)
        target_link_libraries(driver PUBLIC alpaqa::casadi-loader casadi)
        target_compile_definitions(driver PUBLIC ALPAQA_HAVE_CASADI)
    endif()
    if (TARGET alpaqa::cutest-interface)
        target_link_libraries(driver PUBLIC alpaqa::cutest-interface)
        target_compile_definitions(driver PUBLIC ALPAQA_HAVE_CUTEST)
    endif()
    if (TARGET alpaqa::casadi-ocp-loader)
        target_link_libraries(driver PUBLIC alpaqa::casadi-ocp-loader)
    endif()
    if (TARGET alpaqa::ipopt-adapter)
        target_link_libraries(driver PUBLIC alpaqa::ipopt-adapter)
        target_compile_definitions(driver PUBLIC WITH_IPOPT)
    endif()
    if (TARGET alpaqa::lbfgsb-adapter)
        target_link_libraries(driver PUBLIC alpaqa::lbfgsb-adapter)
        target_compile_definitions(driver PUBLIC WITH_LBFGSB)
    endif()
endif()

# Installation
include(cmake/Install.cmake)
