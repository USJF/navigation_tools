#ifndef NAVIGATION_CALCULATOR_HPP
#define NAVIGATION_CALCULATOR_HPP

#include "coordinate_converter.hpp"
#include <vector>
#include <cmath>
#include <iostream>

class NavigationCalculator : public CoordinateConverter
{
    public:
    NavigationCalculator();
    ~NavigationCalculator();

    struct PolarData
    {
        double range;
        double bearing;
    };

    PolarData calculatePolarData(LLA start_point, LLA target_point);

    double calculateRn(double latitude);
    double calculateRe(double latitude);

    private:
    double Rp;
    double R0;
    double e;
    double e2;

    double radiansToDegrees(double radians);
    double degreesToRadians(double degrees);
    double calculateBearing(double delta_east, double delta_north);
};

#endif