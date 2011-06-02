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

