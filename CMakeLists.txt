cmake_minimum_required(VERSION 4.0.2)

project(allignment)

set(CMAKE_CXX_STANDARD 11)

# Add include path for bioparser
include_directories(${CMAKE_SOURCE_DIR}/bioparser/include)

# Add ZLIB so that bioparser works
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

# Define the executable before linking libraries
add_executable(${PROJECT_NAME} main.cpp)

# Now link libraries to the target
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

