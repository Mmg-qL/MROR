#一、描述性命令
cmake_minimum_required(VERSION 2.8.3)
project(MROR_filter)

#二、指定C++标准
SET(CMAKE_CXX_STANDARD 17)

#三、include头文件
include_directories(
			${catkin_INCLUDE_DIRS}
			${PCL_LIBRARIES}
		    ${EIGEN3_INCLUDE_DIR}
			include)

#四、查找库
find_package(PCL 1.8 REQUIRED)
find_package(catkin REQUIRED COMPONENTS
			roscpp
			rospy
			sensor_msgs
			std_msgs
			pcl_conversions)
			
catkin_package(INCLUDE_DIRS include)

#五、添加可执行文件
add_executable(dynamicMROR src/dynamicMROR.cpp)


#六、链接库
add_library(${PROJECT_NAME}_core src/MROR.cpp)
add_dependencies(${PROJECT_NAME}_core ${catkin_EXPORTED_TARGETS})			 

target_link_libraries(${PROJECT_NAME}_core 
					  ${PCL_LIBRARIES}
					  ${catkin_LIBRARIES} 
					  ${EIGEN3_INCLUDE_DIR})

target_link_libraries(dynamicMROR ${catkin_LIBRARIES}
					  ${PCL_LIBRARIES} ${PROJECT_NAME}_core)
					  
add_definitions(${PCL_DEFINITIONS})	
