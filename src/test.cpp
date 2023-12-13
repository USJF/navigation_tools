#include "navigation_tools/navigation_calculator.hpp"
#include <iostream>
int main()
{
    double min_lat = 32.594360;
    double min_lon = -85.296170;
    double max_lat = 32.596404;
    double max_lon = -85.294538;

    double alt = 213;

    double x = 441171.57858859503;
    double y = -5360799.309623537;
    double z = 3416339.7421602733;

    NavigationCalculator::LLA ref;
    ref.latitude = min_lat;
    ref.longitude = min_lon;
    ref.alt = alt;

    NavigationCalculator::ECEF point;
    point.x = x;
    point.y = y;
    point.z = z;

    NavigationCalculator navigation_calculator;
    NavigationCalculator::Vector3 enu = navigation_calculator.convertECEF2ENU(point, ref);
    std::cout<<enu.x<<std::endl;
    std::cout<<enu.y<<std::endl;
    std::cout<<enu.z<<std::endl;
}