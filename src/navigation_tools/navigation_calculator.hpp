#ifndef NAVIGATION_CALCULATOR_HPP
#define NAVIGATION_CALCULATOR_HPP

#include "coordinate_converter.hpp"
#include <vector>
#include <cmath>

class NavigationCalculator : public CoordinateConverter
{
    public:
    NavigationCalculator();
    ~NavigationCalculator();

    struct PolarData
    {
        double displacement;
        double haversine;
        double bearing;
    };

    PolarData calculatePolarData(LLA start_point, LLA target_point);

    double calculateRn(double latitude);
    double calculateRe(double latitude);
    double haversine(double latitude1, double latitude2, double longitude1, double longitude2);

    private:
    double Rp;
    double R0;
    double e;
    double e2;

    double delta_x;
    double delta_y;
    double delta_z;

    // Polar Calculation Variables 
    PolarData output_data;
    double local_rotation_coefs[3][3];
    ECEF ecef_start;
    ECEF ecef_target;
    double latitude;
    double longitude;
    double magnitude;
    double enu[3];

    double radiansToDegrees(double radians);
    double degreesToRadians(double degrees);
    double calculateBearing(double delta_east, double delta_north);
};

#endif