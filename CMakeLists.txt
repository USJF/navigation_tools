cmake_minimum_required(VERSION 3.5)

project(navigation_tools)

include_directories(include)
include_directories(src/navigation_tools/)

set(NAV_SRCS
    src/navigation_tools/navigation_calculator.cpp 
    src/navigation_tools/coordinate_converter.cpp   
)

add_executable(main src/test.cpp ${NAV_SRCS})

