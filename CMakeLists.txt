#一、描述性命令
cmake_minimum_required(VERSION 2.8.3)

#二、指定C++标准
SET(CMAKE_CXX_STANDARD 17)

#三、include头文件
include_directories(${PCL_LIBRARIES}
		    ${EIGEN3_INCLUDE_DIR})

#四、查找库
find_package(PCL 1.8 REQUIRED)

#五、添加可执行文件
add_executable(MROR MROR.cpp)

#六、链接库
target_link_libraries(MROR ${PCL_LIBRARIES} ${EIGEN3_INCLUDE_DIR})
