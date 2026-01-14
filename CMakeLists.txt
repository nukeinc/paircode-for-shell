# cmake_minimum_required(VERSION 3.26)

# # --- 1. 项目名称（自动使用当前文件夹名） ---
# get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
# string(REPLACE " " "_" ProjectId ${ProjectId})
# project(${ProjectId} CXX)

# # --- 2. C++ 标准设置 (C++17) ---
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)



# # --- 3. 输出目录 (编译出的程序会放在 bin 文件夹) ---
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# # --- 4. 查找 OpenMP (Mac/Linux 通用) ---
# find_package(OpenMP REQUIRED)

# # --- 5. 包含头文件目录 ---
# include_directories(inc src)

# # --- 6. 自动扫描并编译所有 .cpp 文件 ---
# add_compile_options(-g -O2 -pg)
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# foreach(file ${files})
#     # 获取文件名（去掉后缀），作为程序名
#     get_filename_component(name ${file} NAME_WE)
    
#     # 添加可执行文件
#     add_executable(${name} ${file})

#     # 链接 OpenMP
#     if(OpenMP_CXX_FOUND)
#         target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX)
#     endif()

#     # 【兼容性修正】Linux 特有配置
#     if(UNIX AND NOT APPLE)
#         target_link_libraries(${name} PUBLIC m)
#     endif()
# endforeach()
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     add_compile_options(-g -O0) # 调试模式：无优化，带符号
# else()
#     add_compile_options(-O3 -g) # 发布模式：全优化
# endif()
# # --- 7. 编译优化选项 (Release 模式生效) ---
# set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -flto")

cmake_minimum_required(VERSION 3.26)

# --- 1. 项目名称 ---
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId} CXX)

# --- 2. C++ 标准 ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- 3. 输出目录 ---
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# --- 4. 查找 OpenMP ---
find_package(OpenMP REQUIRED)

# --- 5. 头文件 ---
include_directories(inc src)

# --- [修改部分] 6. 区分 Debug 和 Release 模式 ---
# 默认设置为 Release 模式，如果你没在命令行指定的话
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Build Type: Debug (含调试符号, 无优化, 开启 gprof)")
    # 调试模式：开启 gprof (-pg)，关闭优化 (-O0)，开启调试信息 (-g)
    add_compile_options(-g -O0 -pg)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
else()
    message(STATUS "Build Type: Release (全速优化, 无 gprof)")
    # 发布模式：开启最高优化 (-O3)，开启架构特定优化 (-march=native)
    # 注意：这里千万不要加 -pg
    add_compile_options(-O3 -march=native -DNDEBUG)
    # 如果你的编译器支持 LTO (Link Time Optimization)，可以取消下面这行的注释
    # add_compile_options(-flto)
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
endif()

# --- 7. 扫描源文件 ---
file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(file ${files})
    get_filename_component(name ${file} NAME_WE)
    add_executable(${name} ${file})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX)
    endif()
    if(TARGET TBB::tbb)
    target_link_libraries(${name} PUBLIC TBB::tbb)
    else()
        # 如果没找到 target，尝试直接链接库名（适配某些旧环境）
        target_link_libraries(${name} PUBLIC tbb)
    endif()

    if(UNIX AND NOT APPLE)
        target_link_libraries(${name} PUBLIC m)
    endif()
endforeach()

# cmake_minimum_required(VERSION 3.26)

# # --- 1. 项目名称 ---
# get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
# string(REPLACE " " "_" ProjectId ${ProjectId})
# project(${ProjectId} CXX)

# # --- 2. C++ 标准 ---
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # --- 3. 输出目录 ---
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# # --- 4. 查找依赖库 ---
# find_package(OpenMP REQUIRED)
# # 查找 MKL (Intel oneAPI 提供的配置)
# find_package(MKL CONFIG REQUIRED)

# # --- 5. 配置 Eigen 使用 MKL 的宏 (关键步骤) ---
# # 定义这个宏后，Eigen 的矩阵乘法等操作会自动调用 MKL
# add_definitions(-DEIGEN_USE_MKL_ALL)

# # --- 6. 头文件 ---
# # 包含 Eigen 和项目本身的头文件
# include_directories(inc src ${EIGEN3_INCLUDE_DIRS})

# # --- 7. 区分 Debug 和 Release 模式 ---
# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE Release)
# endif()

# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     message(STATUS "Build Type: Debug")
#     add_compile_options(-g -O0 -pg)
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# else()
#     message(STATUS "Build Type: Release")
#     # MKL 配合 -O3 和 -march=native 能发挥最大威力
#     add_compile_options(-O3 -march=native -DNDEBUG)
# endif()

# # --- 8. 扫描源文件并链接 ---
# file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# foreach(file ${files})
#     get_filename_component(name ${file} NAME_WE)
#     add_executable(${name} ${file})

#     # 链接 OpenMP
#     target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX)

#     # 链接 MKL (MKL::MKL 会自动处理 TBB 或 OpenMP 的并行链接)
#     target_link_libraries(${name} PUBLIC MKL::MKL)

#     # 链接数学库
#     if(UNIX AND NOT APPLE)
#         target_link_libraries(${name} PUBLIC m)
#     endif()
# endforeach()