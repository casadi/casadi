# my examples
find_package(IPOPT QUIET)

IF(IPOPT_FOUND)
include_directories(${IPOPT_INCLUDE_DIR})

include_directories(../cppocp)

add_executable(gliderCppocp glider.cpp)
target_link_libraries(gliderCppocp
  cppocp
  ipopt_interface 
  casadi 
  ${IPOPT_LIBRARIES} 
)

add_executable(rocketCppocp rocket.cpp)
target_link_libraries(rocketCppocp
  cppocp
  ipopt_interface 
  casadi 
  ${IPOPT_LIBRARIES} 
)

ENDIF(IPOPT_FOUND)



find_package(SNOPT QUIET)

IF(SNOPT_FOUND)
include_directories(../snopt)
include_directories(${SNOPT_INCLUDE_DIRS})

add_executable(rocketSnopt rocketSnopt.cpp)
add_executable(my_sntoya my_sntoya.cpp)
target_link_libraries(rocketSnopt
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
  ${IPOPT_LIBRARIES}
  ipopt_interface
)
target_link_libraries(my_sntoya
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
  ${IPOPT_LIBRARIES}
  ipopt_interface
)


ENDIF(SNOPT_FOUND)
