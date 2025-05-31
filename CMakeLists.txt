cmake_minimum_required(VERSION 4.0.2)

project(allignment)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Add the alignment library
add_library(team_alignment 
    team_alignment/team_alignment.cpp
    team_alignment/team_alignment.hpp)
# Add the minimizers library
add_library(team_minimizers
    team_minimizers/team_minimizers.cpp
    team_minimizers/team_minimizers.hpp)


target_include_directories(team_alignment PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/team_alignment )
target_include_directories(team_minimizers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/team_minimizers)

# Add include path for bioparser
include_directories(${CMAKE_SOURCE_DIR}/bioparser/include)

# Add ZLIB so that bioparser works
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

# Define the executable before linking libraries
add_executable(${PROJECT_NAME} team_mapper.cpp)


# Link the alignment library to the mapper
target_link_libraries(${PROJECT_NAME} PRIVATE team_alignment)
# Link the minimizers library to the mapper
target_link_libraries(${PROJECT_NAME} PRIVATE team_minimizers)
# Now link libraries to the target
target_link_libraries(${PROJECT_NAME} PRIVATE ${ZLIB_LIBRARIES})





