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
target_link_libraries(rocketSnopt
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_executable(gliderSnopt gliderSnopt.cpp)
target_link_libraries(gliderSnopt
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_executable(cartpoleSnopt cartpoleSnopt.cpp)
target_link_libraries(cartpoleSnopt
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_library(rocketmex STATIC
    rocketMexLib.cpp
)
add_executable(raptors raptors.cpp)
target_link_libraries(raptors
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_executable(lqrCartpole lqrCartpole.cpp)
target_link_libraries(lqrCartpole
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_executable(lqrSpring lqrSpring.cpp)
target_link_libraries(lqrSpring
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_executable(ddpDoubleCartpole ddpDoubleCartpole.cpp)
target_link_libraries(ddpDoubleCartpole
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_executable(ddpCartpole ddpCartpole.cpp)
target_link_libraries(ddpCartpole
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

add_executable(ddpSpring ddpSpring.cpp)
target_link_libraries(ddpSpring
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)

ENDIF(SNOPT_FOUND)

find_package(MATLAB)

IF(MATLAB_FOUND)

include_directories(../cppocp)
include_directories(../matlab)
include_directories(${MATLAB_INCLUDE_DIRS})

add_executable(testmatlab testmatlab.cpp)
target_link_libraries(testmatlab
	matlab
)


ENDIF(MATLAB_FOUND)
