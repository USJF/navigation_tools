#include <iostream>

#include "navigation_calculator.hpp"

int main(int argc, char **argv)
{
    NavigationCalculator calculator;

    NavigationCalculator::LLA start_point;
    start_point.latitude = 34.862976;
    start_point.longitude = -86.778458;
    start_point.alt = 135;

    NavigationCalculator::LLA target_point;
    target_point.latitude = 35.042636;
    target_point.longitude = -86.983825;
    target_point.alt = 202;

    NavigationCalculator::PolarData data = calculator.calculatePolarData(start_point, target_point);
    std::cout<<data.range<<"\t"<<data.bearing<<std::endl;
    
}