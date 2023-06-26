#include <iostream>
#include <iomanip>

#include "navigation_calculator.hpp"

int main(int argc, char **argv)
{

    NavigationCalculator calculator;

    NavigationCalculator::LLA point;
    NavigationCalculator::LLA start_point = point;

    point.latitude = 45;
    point.longitude = -120;
    point.alt = 10000;

    double total_distance = 885 * 8 * 1000; // km/h * h * m/km = meters
    //total_distance = 100000000;
    double current_distance = 0;
    double heading = 45 * M_PI / 180;

    std::cout << total_distance << std::endl;
    std::cout << heading << std::endl;

    NavigationCalculator::LLA end_point = calculator.findLocationFromPolarData(point, total_distance, heading);

    while (current_distance < total_distance)
    {
        point = calculator.findLocationFromPolarData(point, 5.0, heading);
        current_distance++;
    }

    std::cout << point.latitude << '\t' << point.longitude << '\t' << point.alt << std::endl;

    std::cout << end_point.latitude << '\t' << end_point.longitude << '\t' << end_point.alt << std::endl;

    std::cout << calculator.haversine(start_point, end_point) << std::endl;   
}