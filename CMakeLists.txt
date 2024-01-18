cmake_minimum_required(VERSION 3.5)

project(navigation_tools)

include_directories(src/navigation_tools/)
add_subdirectory(src/navigation_tools)

add_executable(main src/test.cpp)
target_link_libraries(main NavigationTools)

