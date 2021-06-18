# CMake 最低版本号要求
cmake_minimum_required(VERSION 3.0)
 
# 项目名称
project(algorithm)


# 添加头文件路径
# include_directories("./X-CUBE-AI/App/")
# include_directories("./X-CUBE-AI/Inc/")

# 添加库文件路径
# link_directories("X-CUBE-AI/Lib")

# 设置AI代码文件的路径
# set(AI_FILE "./X-CUBE-AI/App/network_data.c" "./X-CUBE-AI/App/network.c")

# 查找指定目录下的所有.cpp与.h文件 并存放到指定变量名SC_FILES中
FILE(GLOB SC_FILES "*.cpp" "*.h")

# source_group(ai FILES ${AI_FILE})

# 指定生成目标
add_executable(${PROJECT_NAME} ${SC_FILES})

# target_link_libraries(${PROJECT_NAME} NetworkRuntime600_CM7_Keil.lib)