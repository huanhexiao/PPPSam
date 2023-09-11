cmake_minimum_required(VERSION 3.0)

project(PPP-SAM)

set( CMAKE_BUILD_TYPE "Debug" )
set(CMAKE_CXX_FLAGS "-std=c++17")

# yaml-cpp
find_package(yaml-cpp REQUIRED)
#GTSAM
# find_package(GTSAM REQUIRED) # Uses installed package
set(GTSAM_LIBRARIES gtsam)   # TODO: automatic search libs
set(GTSAM_UNSTABLE_LIBRARIES gtsam_unstable)

######  include and library directories.
include_directories(
	${GTSAM_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/PPPSamlib
	${PROJECT_SOURCE_DIR}/rtklibfunc/include
	)

file(GLOB SRC_FILES "rtklibfunc/src/*.cpp" "rtklibfunc/src/*.c" )
add_library(RTKLib4PPPSAM ${SRC_FILES})
add_library(${PROJECT_NAME} SHARED
    ${SRC_FILES}
	${PROJECT_SOURCE_DIR}/PPPSamlib/PPPSamlib.cpp
)

target_link_libraries(${PROJECT_NAME}
    ${YAML_CPP_LIBRARIES}
	${GTSAM_LIBRARIES}
	${GTSAM_UNSTABLE_LIBRARIES}
	)

#############################
# add_executable(pppBayesTree_isam2_rl pppBayesTree_isam2_rl.cpp)
# target_link_libraries(pppBayesTree_isam2_rl ${PROJECT_NAME} RTKLib4PPPSAM)

# add_executable(pppBayesTree_isam2_rl_newstate pppBayesTree_isam2_rl_newstate.cpp)
# target_link_libraries(pppBayesTree_isam2_rl_newstate ${PROJECT_NAME} RTKLib4PPPSAM)

add_executable(PPPSAM PPPSAM.cpp)
target_link_libraries(PPPSAM ${PROJECT_NAME})

add_executable(PPPSAM-dyna PPPSAM-dyna.cpp)
target_link_libraries(PPPSAM-dyna ${PROJECT_NAME})