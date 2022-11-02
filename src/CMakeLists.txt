find_package(Eigen3 REQUIRED)

add_library(alpaqa
    "src/util/type-erasure.cpp"
    "src/util/demangled-typename.cpp"
    "src/util/print.cpp"
    "src/accelerators/lbfgs.cpp"
    "src/problem/problem-counters.cpp"
    "src/outer/alm.cpp"
    "src/inner/panoc.cpp"
    "src/inner/structured-panoc.cpp"
    "src/panoc-alm.cpp"
    "src/structured-panoc-alm.cpp"
    "src/inner/internal/solverstatus.cpp"
    "include/alpaqa/util/type-erasure.hpp"
)
target_compile_features(alpaqa PUBLIC cxx_std_20)
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_definitions(alpaqa PUBLIC
    $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
    $<$<CONFIG:Debug>:EIGEN_RUNTIME_NO_MALLOC>)
target_link_libraries(alpaqa PUBLIC Eigen3::Eigen)
target_link_libraries(alpaqa PRIVATE warnings)
set_target_properties(alpaqa PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                        VISIBILITY_INLINES_HIDDEN true)
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_options(alpaqa PRIVATE "LINKER:--exclude-libs,ALL")
endif()
target_compile_definitions(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:ALPAQA_WITH_QUAD_PRECISION>)
target_link_libraries(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:quadmath>)
add_library(alpaqa::alpaqa ALIAS alpaqa)
list(APPEND ALPAQA_INSTALL_TARGETS alpaqa)

# DLL import/export
include(GenerateExportHeader)
generate_export_header(alpaqa
    EXPORT_FILE_NAME alpaqa/export.h)
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

# CasADi
if (ALPAQA_WITH_CASADI)
    find_package(casadi REQUIRED)
    add_library(casadi-loader 
        "interop/casadi/CasADiLoader.cpp"
        "include/alpaqa/interop/casadi/CasADiFunctionWrapper.hpp"
        "include/alpaqa/interop/casadi/CasADiLoader.hpp")
    find_package(Threads REQUIRED)
    target_link_libraries(casadi-loader 
        PUBLIC
            alpaqa::alpaqa
        PRIVATE
            casadi
            Threads::Threads)
    target_link_libraries(casadi-loader PRIVATE warnings)
    set_target_properties(casadi-loader PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                        VISIBILITY_INLINES_HIDDEN true)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        target_link_options(casadi-loader PRIVATE "LINKER:--exclude-libs,ALL")
    endif()
    add_library(alpaqa::casadi-loader ALIAS casadi-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS casadi-loader)
    # DLL import/export
    generate_export_header(casadi-loader
    EXPORT_FILE_NAME alpaqa/casadi-loader-export.h)
    target_include_directories(casadi-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

    message(STATUS "Compiling WITH CasADi support")
else()
    message(STATUS "Compiling WITHOUT CasADi support")
endif()

# Installation
include(cmake/Install.cmake)
