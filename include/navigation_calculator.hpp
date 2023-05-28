#ifndef NAVIGATION_CALCULATOR_HPP
#define NAVIGATION_CALCULATOR_HPP

#include "coordinate_converter.hpp"

class NavigationCalculator : public CoordinateConverter
{
    public:
    NavigationCalculator();
    ~NavigationCalculator();


    private:
    double Rp;
    double R0;
    double e;
    double e2;

    CoordinateConverter::UTM utm_coords;
    CoordinateConverter::MGRS mgrs_coords;
    CoordinateConverter::LLA lla_coords;
    CoordinateConverter::LLADMS lladms_coords;
    CoordinateConverter::ECEF ecef_coords;
};

#endif