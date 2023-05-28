#include <iostream>

#include "navigation_calculator.hpp"

int main(int argc, char **argv)
{
    /*CoordinateConverter converter;
    CoordinateConverter::LLA lla_coords;
    lla_coords.latitude = 34.862976;
    lla_coords.longitude = -86.778458;
    lla_coords.alt = 186;
    CoordinateConverter::UTM utm_coords = converter.convertLLA2UTM(lla_coords);
    std::cout<<utm_coords.easting<<"\t"<<utm_coords.northing<<"\t"<<utm_coords.grid_zone<<std::endl;
    CoordinateConverter::MGRS mgrs_coords = converter.convertLLA2MGRS(lla_coords);
    std::cout<<mgrs_coords.easting<<"\t"<<mgrs_coords.northing<<"\t"<<mgrs_coords.grid_zone<<mgrs_coords.grid_letter<<mgrs_coords.false_easting<<mgrs_coords.false_northing<<std::endl;*/

    NavigationCalculator calculator;

    NavigationCalculator::LLA lla_coords;
    lla_coords.latitude = 34.862976;
    lla_coords.longitude = -86.778458;
    lla_coords.alt = 186;

    NavigationCalculator::UTM utm_coords = calculator.convertLLA2UTM(lla_coords);
    std::cout<<utm_coords.easting<<"\t"<<utm_coords.northing<<"\t"<<utm_coords.grid_zone<<std::endl;
    NavigationCalculator::MGRS mgrs_coords = calculator.convertLLA2MGRS(lla_coords);
    std::cout<<mgrs_coords.easting<<"\t"<<mgrs_coords.northing<<"\t"<<mgrs_coords.grid_zone<<mgrs_coords.grid_letter<<mgrs_coords.false_easting<<mgrs_coords.false_northing<<std::endl;
}