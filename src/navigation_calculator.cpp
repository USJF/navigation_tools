#include "navigation_calculator.hpp"

NavigationCalculator::NavigationCalculator()
{
    // Earth Modeling Constants
    R0 = 6378137.0;
    Rp = 6356752.3142;
    e2 = 1 - ((Rp*Rp)/(R0 * R0));
    e = sqrt(e2);
}

NavigationCalculator::~NavigationCalculator()
{

}

double NavigationCalculator::calculateRn(double latitude)
{
    double num = (1 - e2) * R0;

    double den = pow(1 - (e2 * pow(sin(latitude), 2)), 3/2);
    return num/den;
}

double NavigationCalculator::calculateRe(double latitude)
{
    double den = sqrt(1 - (e2 * pow(sin(latitude), 2)));
    return R0 / den;
}

NavigationCalculator::PolarData NavigationCalculator::calculatePolarData(LLA start_point, LLA target_point)
{
    latitude = degreesToRadians(start_point.latitude);
    longitude = degreesToRadians(start_point.longitude);

    // Convert Points to ECEF
    ecef_start = convertLLA2ECEF(start_point);
    ecef_target = convertLLA2ECEF(target_point);

    // Calculate Vector Characteristics in ECEF
    delta_x = ecef_target.x - ecef_start.x;
    delta_y = ecef_target.y - ecef_start.y;
    delta_z = ecef_target.z - ecef_start.z;
    magnitude = sqrt(pow(delta_x, 2) + pow(delta_y, 2) + pow(delta_z, 2));

    // Create ECEF to ENU Rotation Matrix
    local_rotation_coefs[0][0] = -1 * sin(longitude);
    local_rotation_coefs[0][1] = cos(longitude);
    local_rotation_coefs[0][2] = 0;
    local_rotation_coefs[1][0] = -1 * cos(longitude) * sin(latitude);
    local_rotation_coefs[1][1] = -1 * sin(longitude) * sin(latitude);
    local_rotation_coefs[1][2] = cos(latitude);
    local_rotation_coefs[2][0] = cos(longitude) * cos(latitude);
    local_rotation_coefs[2][1] = sin(longitude) * cos(latitude);
    local_rotation_coefs[2][2] = sin(latitude);

    // Rotate ECEF Unit Vector to ENU
    enu[0] = delta_x/magnitude * local_rotation_coefs[0][0] + delta_y/magnitude * local_rotation_coefs[0][1] + delta_z/magnitude * local_rotation_coefs[0][2];
    enu[1] = delta_x/magnitude * local_rotation_coefs[1][0] + delta_y/magnitude * local_rotation_coefs[1][1] + delta_z/magnitude * local_rotation_coefs[1][2];
    enu[2] = delta_x/magnitude * local_rotation_coefs[2][0] + delta_y/magnitude * local_rotation_coefs[2][1] + delta_z/magnitude * local_rotation_coefs[2][2];

    // Store Outputs
    output_data.displacement = magnitude;
    output_data.bearing = calculateBearing(enu[0], enu[1]);

    output_data.haversine = haversine(start_point.latitude, target_point.latitude, start_point.longitude, target_point.longitude);

    if (output_data.bearing < 0)
    {
        output_data.bearing += 360;
    }

    return output_data;
}

double NavigationCalculator::degreesToRadians(double degrees)
{
    return degrees * M_PI / 180.0;
}

double NavigationCalculator::radiansToDegrees(double radians)
{
    return radians * 180.0 / M_PI;
}

double NavigationCalculator::calculateBearing(double delta_east,  double delta_north)
{
    // Calculate ENU Bearing from Rotated Unit Vector
    double target_bearing = atan(delta_east/delta_north);
    if (delta_east<0)
    {
        if (delta_north<0)
        {
            target_bearing += M_PI;
        }
        else
        {
            target_bearing += 2*M_PI;
        }
    }
    else
    {
        if (delta_north<0)
        {
            target_bearing += M_PI;
        }
        else
        {
            //do nothing
        }
    }
    
    double relative_bearing = target_bearing;
    if(relative_bearing>M_PI)
    {
        relative_bearing += -2*M_PI;
    }

    return radiansToDegrees(relative_bearing);
}

double NavigationCalculator::haversine(double latitude1, double latitude2, double longitude1, double longitude2)
{
    // Ground Distance Calculation Using Haversine Formula
    latitude1 = degreesToRadians(latitude1);
    latitude2 = degreesToRadians(latitude2);
    double delta_latitude = latitude2 - latitude1;

    double delta_longitude = degreesToRadians(longitude2 - longitude1);
    double a = pow(sin(delta_latitude/2), 2) + cos(latitude1) * cos(latitude2) * pow(sin(delta_longitude/2), 2);

    return R0 * 2 * atan2(sqrt(a), sqrt(1-a));
}