cmake_minimum_required(VERSION 3.10)
project(NuclearPairApproximation VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Library
add_library(pairapprox STATIC
    src/PairApproximation.cpp
)

# Main executable
add_executable(paircode src/main.cpp)
target_link_libraries(paircode pairapprox)

# Testing
enable_testing()

# Test executable
add_executable(test_pair tests/test_pair.cpp)
target_link_libraries(test_pair pairapprox)

add_test(NAME TestPairApproximation COMMAND test_pair)

# Install
install(TARGETS paircode RUNTIME DESTINATION bin)
install(TARGETS pairapprox ARCHIVE DESTINATION lib)
install(FILES include/PairApproximation.hpp DESTINATION include)
