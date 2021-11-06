# Choose problem set
# ------------------

# include(AllProblems.cmake)
# include(BoxProblems.cmake)
include(BoxProblemsNoQP.cmake)

find_package(yaml-cpp)
if (yaml-cpp_FOUND AND CUTEst_FOUND)
    add_library(yaml-encoder src/YAMLEncoder.cpp)
    target_include_directories(yaml-encoder
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
        PRIVATE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    )
    target_link_libraries(yaml-encoder PUBLIC yaml-cpp cutest-loader alpaqa)
    
    string(JOIN "\", \"" CUTEst_PROBLEM_LIST_STRINGS ${CUTEst_PROBLEM_LIST})
    string(CONCAT CUTEst_PROBLEM_LIST_STRINGS "\"" "${CUTEst_PROBLEM_LIST_STRINGS}" "\"")
    configure_file(include/drivers/CUTEstProblemList.h.in
                   include/drivers/CUTEstProblemList.h @ONLY)

    set(TESTNAME "strucpanoc-21-10-hessheur10")
    set(TESTRESULT_DIR ${CMAKE_BINARY_DIR}/testresults/${TESTNAME}/CUTEst)
    file(MAKE_DIRECTORY ${TESTRESULT_DIR})
    foreach(P ${CUTEst_PROBLEM_LIST})
        cutest_sif_problem(${P})
    endforeach()

    add_executable(cutest-driver src/driver.cpp)
    target_compile_definitions(cutest-driver PRIVATE "SOLVER=SOLVER_PANOC_2ND_LBFGS")
    target_link_libraries(cutest-driver PRIVATE yaml-encoder)

    foreach(P ${CUTEst_PROBLEM_LIST})
        add_custom_command(OUTPUT ${TESTRESULT_DIR}/${P}.yaml
                           COMMAND cutest-driver ${P} ${TESTRESULT_DIR}
                           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                           DEPENDS CUTEst::problem-${P} cutest-driver)
        add_custom_target(cutest-result-${TESTNAME}-${P}
                          DEPENDS ${TESTRESULT_DIR}/${P}.yaml)
        list(APPEND CUTEst_PROBLEMS_RESULT_TARGETS 
                    cutest-result-${TESTNAME}-${P})
    endforeach()
    add_custom_command(OUTPUT ${TESTRESULT_DIR}/parameters.yaml
                       COMMAND cutest-driver parameters ${TESTRESULT_DIR}
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS cutest-driver)
    add_custom_target(cutest-result-${TESTNAME}-parameters
                      DEPENDS ${TESTRESULT_DIR}/parameters.yaml)
    list(APPEND CUTEst_PROBLEMS_RESULT_TARGETS 
                cutest-result-${TESTNAME}-parameters)
    message(STATUS cutest-results-${TESTNAME}-export)
    add_custom_target(cutest-results-${TESTNAME}
                      DEPENDS ${CUTEst_PROBLEMS_RESULT_TARGETS})
    cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)
    add_custom_target(cutest-results-${TESTNAME}-export
        COMMAND ${CMAKE_COMMAND} -E copy_directory 
            ${CMAKE_BINARY_DIR}/testresults/${TESTNAME}
            ${CMAKE_CURRENT_LIST_DIR}/../testresults/${HOSTNAME}/${TESTNAME}
        DEPENDS cutest-results-${TESTNAME})
endif()
