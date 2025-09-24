cmake_minimum_required(VERSION 3.26)
# 项目名称（用当前目录名）
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId} CXX)
# C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
# 输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
# 查找 OpenMP（Homebrew LLVM 会找到 libomp）
find_package(OpenMP REQUIRED)
# 包含头文件目录
include_directories(inc src)
# 获取主目录下所有 .cpp 文件
file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
# 为每个 .cpp 文件生成一个可执行文件
foreach(file ${files})
    get_filename_component(name ${file} NAME_WE)   # 去掉扩展名
    add_executable(${name} ${file})
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX)
    endif()
endforeach()
# 编译优化
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto")
#set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")