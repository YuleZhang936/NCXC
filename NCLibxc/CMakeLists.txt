cmake_minimum_required(VERSION 3.10)
project(NCLibxc)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 启用 O3 优化
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# 添加头文件搜索路径
include_directories(/public1/home/scg0213/software-scg0213/libxc-7.0.0/install/include)

# 添加库文件搜索路径
link_directories(/public1/home/scg0213/software-scg0213/libxc-7.0.0/install/lib)

# 定义源文件列表，排除 test_NCLibxc.cpp
set(SOURCE_FILES
    NCLibxc.cpp
    LebedevGrid.cpp
    interface_to_libxc.cpp
    FibonacciGrid.cpp
    lda.cpp
    gga.cpp
    mgga.cpp
    math.cpp
    print.cpp
    benchmark_tests.cpp
    torque.cpp
    gga_col.cpp
)

# 生成动态库
add_library(NCLibxc SHARED ${SOURCE_FILES})

# 链接 libxc 库
target_link_libraries(NCLibxc xc)

# 安装动态库和头文件
install(TARGETS NCLibxc DESTINATION lib)
install(DIRECTORY . DESTINATION include FILES_MATCHING PATTERN "*.h")