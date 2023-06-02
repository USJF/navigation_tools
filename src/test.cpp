#include <iostream>
#include <iomanip>

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

    NavigationCalculator::LLA lla;
    lla.latitude = 32.436966;
    lla.longitude = -94.797127;
    lla.alt = 67.0;
    NavigationCalculator::ECEF ecef = calculator.convertLLA2ECEF(lla);
    lla = calculator.convertECEF2LLA(ecef);
    
    std::cout<<std::fixed<<std::setprecision(6)<<std::endl;
    std::cout<<lla.latitude<<"\t"<<lla.longitude<<"\t"<<lla.alt<<std::endl;

}