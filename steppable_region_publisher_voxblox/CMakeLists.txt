cmake_minimum_required(VERSION 3.0.2)
project(steppable_region_publisher_voxblox)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED COMPONENTS
  voxblox
  roslib
  choreonoid_viewer
  choreonoid_voxblox
  )

find_package(choreonoid REQUIRED)

###################################
## catkin specific configuration ##
###################################
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES steppable_region_publisher_voxblox
  CATKIN_DEPENDS voxblox
#  DEPENDS system_lib
)

###########
## Build ##
###########

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${CHOREONOID_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  src/steppable_region_publisher_voxblox.cpp
  )

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${CHOREONOID_BODY_LIBRARIES}
  )

add_library(CnoidSteppableRegionPublisherVoxBloxPlugin SHARED
  src/CnoidSteppableRegionPublisherVoxBloxPlugin.cpp
  src/sample1_vnoid.cpp
)
target_link_libraries(CnoidSteppableRegionPublisherVoxBloxPlugin
  ${catkin_LIBRARIES}
  ${CHOREONOID_BODY_PLUGIN_LIBRARIES}
  ${PROJECT_NAME}
  )
set_target_properties(CnoidSteppableRegionPublisherVoxBloxPlugin PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CHOREONOID_PLUGIN_DIR})

#############
## Install ##
#############

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h"
)
