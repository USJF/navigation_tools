set(NAV_SRCS
    navigation_calculator.cpp 
    coordinate_converter.cpp   
)

add_library(NavigationTools SHARED ${NAV_SRCS})