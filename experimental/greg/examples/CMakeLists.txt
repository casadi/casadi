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

add_executable(raptors raptors.cpp)
target_link_libraries(raptors
  cppocp
  casadi 
  ${SNOPT_LIBRARIES}
  snopt_interface
)


#add_executable(my_sntoya my_sntoya.cpp)
#target_link_libraries(my_sntoya
#  cppocp
#  casadi 
#  ${SNOPT_LIBRARIES}
#  snopt_interface
#)


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