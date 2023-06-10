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

std::string NavigationCalculator::getLocationData()
{
    output_string = "";
    std::string temp;
    //MGRS
    output_string = std::to_string(location_mgrs.grid_zone);
    output_string.append("/");
    temp = location_mgrs.grid_letter + location_mgrs.false_easting + location_mgrs.false_northing;
    output_string.append(temp);
    output_string.append("/");
    output_string.append(std::to_string(location_mgrs.easting));
    output_string.append("/");
    output_string.append(std::to_string(location_mgrs.northing));
    output_string.append("/");
    // UTM
    temp = std::to_string(location_utm.grid_zone) + "/" + std::to_string(location_utm.easting) + "/" + std::to_string(location_utm.northing) + "/";
    output_string.append(temp);
    //LLA
    temp = std::to_string(location_lla.latitude) + "/" + std::to_string(location_lla.longitude) + "/";
    output_string.append(temp);
    // LLADMS
    temp = std::to_string(location_lladms.latitude_degrees) + "/" + std::to_string(abs(location_lladms.latitude_minutes)) + "/" + std::to_string(fabs(location_lladms.latitude_seconds)) + "/";
    output_string.append(temp);
    temp = std::to_string(location_lladms.longitude_degrees) + "/" + std::to_string(abs(location_lladms.longitude_minutes)) + "/" + std::to_string(fabs(location_lladms.longitude_seconds)) + "/";
    output_string.append(temp);
    //ECEF
    temp = std::to_string(location_ecef.x) + "/" + std::to_string(location_ecef.y) + "/" + std::to_string(location_ecef.z);
    output_string.append(temp);

    return output_string;
}

// Set Location to ECEF Coodinates
void NavigationCalculator::setLocationECEF(double x, double y, double z)
{
    location_ecef.x = x;
    location_ecef.y = y;
    location_ecef.z = z;

    location_lla = convertECEF2LLA(location_ecef);
    location_utm = convertLLA2UTM(location_lla);
    location_mgrs = convertUTM2MGRS(location_utm);
    location_lladms = convertLLA2LLADMS(location_lla);
}

void NavigationCalculator::setLocationECEF(ECEF ecef_coordinates)
{
    location_ecef = ecef_coordinates;
    location_lla = convertECEF2LLA(location_ecef);
    location_utm = convertLLA2UTM(location_lla);
    location_mgrs = convertUTM2MGRS(location_utm);
    location_lladms = convertLLA2LLADMS(location_lla);
}

// Set Location to LLA Coordinates
void NavigationCalculator::setLocationLLA(double latitude, double longitude, double alt)
{
    location_lla.latitude = latitude;
    location_lla.longitude = longitude;
    location_lla.alt = alt;

    location_ecef = convertLLA2ECEF(location_lla);
    location_lladms = convertLLA2LLADMS(location_lla);
    location_utm = convertLLA2UTM(location_lla);
    location_mgrs = convertUTM2MGRS(location_utm);
}

void NavigationCalculator::setLocationLLA(LLA lla_coordinates)
{
    location_lla = lla_coordinates;
    location_ecef = convertLLA2ECEF(location_lla);
    location_lladms = convertLLA2LLADMS(location_lla);
    location_utm = convertLLA2UTM(location_lla);
    location_mgrs = convertUTM2MGRS(location_utm);
}

// Set Location to LLADMS Coordinates
void NavigationCalculator::setLocationLLADMS(int latitude_degrees, int latitude_minutes, double latitude_seconds, int longitude_degrees, int longitude_minutes, double longitude_seconds, double alt)
{
    location_lladms.latitude_degrees = latitude_degrees;
    location_lladms.latitude_minutes = latitude_minutes;
    location_lladms.latitude_seconds = latitude_seconds;
    location_lladms.longitude_degrees = longitude_degrees;
    location_lladms.longitude_minutes = longitude_minutes;
    location_lladms.longitude_seconds = longitude_seconds;
    location_lladms.alt = alt;

    location_lla = convertLLADMS2LLA(location_lladms);
    location_ecef = convertLLA2ECEF(location_lla);
    location_utm = convertLLA2UTM(location_lla);
    location_mgrs = convertUTM2MGRS(location_utm);
}

void NavigationCalculator::setLocationLLADMS(LLADMS lladms_coordinates)
{
    location_lladms = lladms_coordinates;
    location_lla = convertLLADMS2LLA(location_lladms);
    location_ecef = convertLLA2ECEF(location_lla);
    location_utm = convertLLA2UTM(location_lla);
    location_mgrs = convertUTM2MGRS(location_utm);
}

// Set Location to MGRS Coordinates
void NavigationCalculator::setLocationMGRS(double easting, double northing, double alt, int grid_number, std::string grid_letter, std::string false_easting, std::string false_northing)
{
    location_mgrs.easting = easting;
    location_mgrs.northing = northing;
    location_mgrs.alt = alt;
    location_mgrs.grid_zone = grid_number;
    location_mgrs.grid_letter = grid_letter;
    location_mgrs.false_easting = false_easting;
    location_mgrs.false_northing = false_northing;

    location_utm = convertMGRS2UTM(location_mgrs);
    location_lla = convertUTM2LLA(location_utm);
    location_lladms = convertLLA2LLADMS(location_lla);
    location_ecef = convertLLA2ECEF(location_lla);
}

void NavigationCalculator::setLocationMGRS(double easting, double northing, double alt, int grid_number, std::string zone_letters)
{
    location_mgrs.easting = easting;
    location_mgrs.northing = northing;
    location_mgrs.alt = alt;
    location_mgrs.grid_zone = grid_number;
    location_mgrs.grid_letter = zone_letters[0];
    location_mgrs.false_easting = zone_letters[1];
    location_mgrs.false_northing = zone_letters[2];

    location_utm = convertMGRS2UTM(location_mgrs);
    location_lla = convertUTM2LLA(location_utm);
    location_lladms = convertLLA2LLADMS(location_lla);
    location_ecef = convertLLA2ECEF(location_lla);
}

void NavigationCalculator::setLocationMGRS(MGRS mgrs_coordinates)
{
    location_mgrs = mgrs_coordinates;
    location_utm = convertMGRS2UTM(location_mgrs);
    location_lla = convertUTM2LLA(location_utm);
    location_lladms = convertLLA2LLADMS(location_lla);
    location_ecef = convertLLA2ECEF(location_lla);
}

// Set Location to UTM Coordinates
void NavigationCalculator::setLocationUTM(double easting, double northing, double alt, int grid_zone)
{
    location_utm.easting = easting;
    location_utm.northing = northing;
    location_utm.alt = alt;
    location_utm.grid_zone = grid_zone;

    location_mgrs = convertUTM2MGRS(location_utm);
    location_lla = convertUTM2LLA(location_utm);
    location_lladms = convertLLA2LLADMS(location_lla);
    location_ecef = convertLLA2ECEF(location_lla);
}

void NavigationCalculator::setLocationUTM(UTM utm_coordinates)
{
    location_utm = utm_coordinates;
    location_mgrs = convertUTM2MGRS(location_utm);
    location_lla = convertUTM2LLA(location_utm);
    location_lladms = convertLLA2LLADMS(location_lla);
    location_ecef = convertLLA2ECEF(location_lla);
}

