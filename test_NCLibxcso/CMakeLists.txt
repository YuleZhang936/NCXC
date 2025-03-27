cmake_minimum_required(VERSION 3.10)
project(TestNCLibxc)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 不启用优化
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

# 添加头文件搜索路径
include_directories(/usr/local/include)

# 添加库文件搜索路径
link_directories(/usr/local/lib)

# 生成可执行文件
add_executable(test_program main.cpp)

# 链接 NCLibxc 动态库
target_link_libraries(test_program NCLibxc)