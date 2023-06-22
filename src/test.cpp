#include <iostream>
#include <iomanip>

#include "navigation_calculator.hpp"

int main(int argc, char **argv)
{

    NavigationCalculator calculator;

    NavigationCalculator::LLA point;

    point.latitude = 45;
    point.longitude = -120;
    point.alt = 10000;

    double total_distance = 885 * 8 * 1000; // km/h * h * m/km = meters
    //total_distance = 100000000;
    double current_distance = 0;
    double heading = 45 * M_PI / 180;

    std::cout << total_distance << std::endl;
    std::cout << heading << std::endl;

    while (current_distance < total_distance)
    {
        point = calculator.findLocationFromPolarData(point, 1.0, heading);
        current_distance++;
    }

    std::cout << point.latitude << '\n' << point.longitude << '\n' << point.alt << std::endl;

}