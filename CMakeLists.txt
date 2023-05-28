cmake_minimum_required(VERSION 3.5)

project(navigation_tools)

include_directories(include)

add_executable(main src/test.cpp src/coordinate_converter.cpp src/navigation_calculator.cpp)

