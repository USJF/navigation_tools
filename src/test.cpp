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

   NavigationCalculator::Vector3 ecef_vector;
   ecef_vector.x = -1;
   ecef_vector.y = 1;
   ecef_vector.z = 0;
   NavigationCalculator::LLA reference;
   reference.latitude = 0;
   reference.longitude = 45;
   reference.alt = 0;
   NavigationCalculator::Vector3 local_vector = calculator.rotateECEFVector2ENU(reference, ecef_vector);
   std::cout<<local_vector.x << "\t"<< local_vector.y <<'\t' << local_vector.z<<std::endl;
}