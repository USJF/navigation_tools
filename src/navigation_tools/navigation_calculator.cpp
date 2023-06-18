#include "navigation_calculator.hpp"

NavigationCalculator::NavigationCalculator()
{
    // Earth Modeling Constants
    R0 = 6378137.0;
    Rp = 6356752.3142;
    e2 = 1 - ((Rp * Rp) / (R0 * R0));
    e = sqrt(e2);

    setLocationLLA(0.0, 0.0, 0.0);
}

NavigationCalculator::NavigationCalculator(LLA lla_initial_position)
{
    // Earth Modeling Constants
    R0 = 6378137.0;
    Rp = 6356752.3142;
    e2 = 1 - ((Rp * Rp) / (R0 * R0));
    e = sqrt(e2);

    setLocationLLA(lla_initial_position);
}

NavigationCalculator::NavigationCalculator(ECEF ecef_initial_position)
{
    // Earth Modeling Constants
    R0 = 6378137.0;
    Rp = 6356752.3142;
    e2 = 1 - ((Rp * Rp) / (R0 * R0));
    e = sqrt(e2);

    setLocationECEF(ecef_initial_position);
}

NavigationCalculator::NavigationCalculator(UTM utm_initial_position)
{
    // Earth Modeling Constants
    R0 = 6378137.0;
    Rp = 6356752.3142;
    e2 = 1 - ((Rp * Rp) / (R0 * R0));
    e = sqrt(e2);

    setLocationUTM(utm_initial_position);
}

NavigationCalculator::NavigationCalculator(LLADMS lladms_initial_position)
{
    // Earth Modeling Constants
    R0 = 6378137.0;
    Rp = 6356752.3142;
    e2 = 1 - ((Rp * Rp) / (R0 * R0));
    e = sqrt(e2);

    setLocationLLADMS(lladms_initial_position);
}

NavigationCalculator::NavigationCalculator(MGRS mgrs_initial_position)
{
    // Earth Modeling Constants
    R0 = 6378137.0;
    Rp = 6356752.3142;
    e2 = 1 - ((Rp * Rp) / (R0 * R0));
    e = sqrt(e2);

    setLocationMGRS(mgrs_initial_position);
}

NavigationCalculator::~NavigationCalculator()
{
}

double NavigationCalculator::calculateRn(double latitude)
{
    double num = (1 - e2) * R0;

    double den = pow(1 - (e2 * pow(sin(latitude), 2)), 3 / 2);
    return num / den;
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
    enu[0] = delta_x / magnitude * local_rotation_coefs[0][0] + delta_y / magnitude * local_rotation_coefs[0][1] + delta_z / magnitude * local_rotation_coefs[0][2];
    enu[1] = delta_x / magnitude * local_rotation_coefs[1][0] + delta_y / magnitude * local_rotation_coefs[1][1] + delta_z / magnitude * local_rotation_coefs[1][2];
    enu[2] = delta_x / magnitude * local_rotation_coefs[2][0] + delta_y / magnitude * local_rotation_coefs[2][1] + delta_z / magnitude * local_rotation_coefs[2][2];

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

double NavigationCalculator::calculateBearing(double delta_east, double delta_north)
{
    // Calculate ENU Bearing from Rotated Unit Vector
    double target_bearing = atan(delta_east / delta_north);
    if (delta_east < 0)
    {
        if (delta_north < 0)
        {
            target_bearing += M_PI;
        }
        else
        {
            target_bearing += 2 * M_PI;
        }
    }
    else
    {
        if (delta_north < 0)
        {
            target_bearing += M_PI;
        }
        else
        {
            // do nothing
        }
    }

    double relative_bearing = target_bearing;
    if (relative_bearing > M_PI)
    {
        relative_bearing += -2 * M_PI;
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
    double a = pow(sin(delta_latitude / 2), 2) + cos(latitude1) * cos(latitude2) * pow(sin(delta_longitude / 2), 2);

    return R0 * 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
}

std::string NavigationCalculator::getLocationData()
{
    output_string = "";
    std::string temp;
    // MGRS
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
    // LLA
    temp = std::to_string(location_lla.latitude) + "/" + std::to_string(location_lla.longitude) + "/";
    output_string.append(temp);
    // LLADMS
    temp = std::to_string(location_lladms.latitude_degrees) + "/" + std::to_string(abs(location_lladms.latitude_minutes)) + "/" + std::to_string(fabs(location_lladms.latitude_seconds)) + "/";
    output_string.append(temp);
    temp = std::to_string(location_lladms.longitude_degrees) + "/" + std::to_string(abs(location_lladms.longitude_minutes)) + "/" + std::to_string(fabs(location_lladms.longitude_seconds)) + "/";
    output_string.append(temp);
    // ECEF
    temp = std::to_string(location_ecef.x) + "/" + std::to_string(location_ecef.y) + "/" + std::to_string(location_ecef.z) + "/" + std::to_string(location_lla.alt);
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

NavigationCalculator::Vector3 NavigationCalculator::calculateLocalRPV(LLA start_point, LLA target_point)
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
    enu[0] = delta_x * local_rotation_coefs[0][0] + delta_y * local_rotation_coefs[0][1] + delta_z * local_rotation_coefs[0][2];
    enu[1] = delta_x * local_rotation_coefs[1][0] + delta_y * local_rotation_coefs[1][1] + delta_z * local_rotation_coefs[1][2];
    enu[2] = delta_x * local_rotation_coefs[2][0] + delta_y * local_rotation_coefs[2][1] + delta_z * local_rotation_coefs[2][2];

    Vector3 local_rpv;
    local_rpv.x = enu[0];
    local_rpv.y = enu[1];
    local_rpv.z = enu[2];

    return local_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::calculateLocalRPV(ECEF start_point, ECEF target_point)
{
    LLA lla_coordinates = convertECEF2LLA(start_point);
    latitude = degreesToRadians(lla_coordinates.latitude);
    longitude = degreesToRadians(lla_coordinates.longitude);

    // Convert Points to ECEF
    ecef_start = start_point;
    ecef_target = target_point;

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
    enu[0] = delta_x * local_rotation_coefs[0][0] + delta_y * local_rotation_coefs[0][1] + delta_z * local_rotation_coefs[0][2];
    enu[1] = delta_x * local_rotation_coefs[1][0] + delta_y * local_rotation_coefs[1][1] + delta_z * local_rotation_coefs[1][2];
    enu[2] = delta_x * local_rotation_coefs[2][0] + delta_y * local_rotation_coefs[2][1] + delta_z * local_rotation_coefs[2][2];

    Vector3 local_rpv;
    local_rpv.x = enu[0];
    local_rpv.y = enu[1];
    local_rpv.z = enu[2];

    return local_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::calculateLocalRPV(LLADMS start_point, LLADMS target_point)
{
    LLA lla_coordinates = convertLLADMS2LLA(start_point);

    latitude = degreesToRadians(lla_coordinates.latitude);
    longitude = degreesToRadians(lla_coordinates.longitude);

    // Convert Points to ECEF
    ecef_start = convertLLA2ECEF(lla_coordinates);
    ecef_target = convertLLA2ECEF(convertLLADMS2LLA(target_point));

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
    enu[0] = delta_x * local_rotation_coefs[0][0] + delta_y * local_rotation_coefs[0][1] + delta_z * local_rotation_coefs[0][2];
    enu[1] = delta_x * local_rotation_coefs[1][0] + delta_y * local_rotation_coefs[1][1] + delta_z * local_rotation_coefs[1][2];
    enu[2] = delta_x * local_rotation_coefs[2][0] + delta_y * local_rotation_coefs[2][1] + delta_z * local_rotation_coefs[2][2];

    Vector3 local_rpv;
    local_rpv.x = enu[0];
    local_rpv.y = enu[1];
    local_rpv.z = enu[2];

    return local_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::calculateLocalRPV(UTM start_point, UTM target_point)
{
    LLA lla_start_point = convertUTM2LLA(start_point);
    LLA lla_target_point = convertUTM2LLA(target_point);

    latitude = degreesToRadians(lla_start_point.latitude);
    longitude = degreesToRadians(lla_start_point.longitude);

    // Convert Points to ECEF
    ecef_start = convertLLA2ECEF(lla_start_point);
    ecef_target = convertLLA2ECEF(lla_target_point);

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
    enu[0] = delta_x * local_rotation_coefs[0][0] + delta_y * local_rotation_coefs[0][1] + delta_z * local_rotation_coefs[0][2];
    enu[1] = delta_x * local_rotation_coefs[1][0] + delta_y * local_rotation_coefs[1][1] + delta_z * local_rotation_coefs[1][2];
    enu[2] = delta_x * local_rotation_coefs[2][0] + delta_y * local_rotation_coefs[2][1] + delta_z * local_rotation_coefs[2][2];

    Vector3 local_rpv;
    local_rpv.x = enu[0];
    local_rpv.y = enu[1];
    local_rpv.z = enu[2];

    return local_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::calculateLocalRPV(MGRS start_point, MGRS target_point)
{
    LLA lla_start_point = convertMGRS2LLA(start_point);
    LLA lla_target_point = convertMGRS2LLA(target_point);

    latitude = degreesToRadians(lla_start_point.latitude);
    longitude = degreesToRadians(lla_start_point.longitude);

    // Convert Points to ECEF
    ecef_start = convertLLA2ECEF(lla_start_point);
    ecef_target = convertLLA2ECEF(lla_target_point);

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
    enu[0] = delta_x * local_rotation_coefs[0][0] + delta_y * local_rotation_coefs[0][1] + delta_z * local_rotation_coefs[0][2];
    enu[1] = delta_x * local_rotation_coefs[1][0] + delta_y * local_rotation_coefs[1][1] + delta_z * local_rotation_coefs[1][2];
    enu[2] = delta_x * local_rotation_coefs[2][0] + delta_y * local_rotation_coefs[2][1] + delta_z * local_rotation_coefs[2][2];

    Vector3 local_rpv;
    local_rpv.x = enu[0];
    local_rpv.y = enu[1];
    local_rpv.z = enu[2];

    return local_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(LLA lla_target_coordinates, double heading)
{
    Vector3 local_rpv = calculateLocalRPV(location_lla, lla_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(LLA lla_target_coordinates, double heading, LLA lla_start_coordinates)
{
    Vector3 local_rpv = calculateLocalRPV(lla_start_coordinates, lla_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(ECEF ecef_target_coordinates, double heading)
{
    Vector3 local_rpv = calculateLocalRPV(location_ecef, ecef_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(ECEF ecef_target_coordinates, double heading, ECEF ecef_start_coordinates)
{
    Vector3 local_rpv = calculateLocalRPV(ecef_start_coordinates, ecef_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(LLADMS lladms_target_coordiantes, double heading)
{
    Vector3 local_rpv = calculateLocalRPV(location_lladms, lladms_target_coordiantes);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(LLADMS lladms_target_coordinates, double heading, LLADMS lladms_start_coordinates)
{
    Vector3 local_rpv = calculateLocalRPV(lladms_start_coordinates, lladms_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(UTM utm_target_coodinates, double heading)
{
    Vector3 local_rpv = calculateLocalRPV(location_utm, utm_target_coodinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(UTM utm_target_coordinates, double heading, UTM utm_start_coordinates)
{
    Vector3 local_rpv = calculateLocalRPV(utm_start_coordinates, utm_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(MGRS mgrs_target_coordinates, double heading)
{
    Vector3 local_rpv = calculateLocalRPV(location_mgrs, mgrs_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::Vector3 NavigationCalculator::getBodyRPV(MGRS mgrs_target_coordinates, double heading, MGRS mgrs_start_coordinates)
{
    Vector3 local_rpv = calculateLocalRPV(mgrs_start_coordinates, mgrs_target_coordinates);
    Vector3 body_rpv;

    body_rpv.x = sin(heading) * local_rpv.x + cos(heading) * local_rpv.y;
    body_rpv.y = cos(heading) * local_rpv.x - sin(heading) * local_rpv.y;
    body_rpv.z = -1 * local_rpv.z;

    return body_rpv;
}

NavigationCalculator::LLA NavigationCalculator::calculateLocation(LLA lla_start_coordinates, double distance, double bearing)
{
    LLA lla_target_coordinates;
    LLA reference_coordinates = lla_start_coordinates;
    double delta;

    lla_start_coordinates.latitude = degreesToRadians(lla_start_coordinates.latitude);
    lla_start_coordinates.longitude = degreesToRadians(lla_start_coordinates.longitude);

    delta = distance / R0;
    lla_target_coordinates.latitude = asin(sin(lla_start_coordinates.latitude) * cos(delta) + cos(lla_start_coordinates.latitude) * sin(delta) * cos(bearing));
    lla_target_coordinates.longitude = lla_start_coordinates.longitude + atan2(sin(bearing) * sin(delta) * cos(lla_start_coordinates.latitude), cos(delta) - sin(lla_start_coordinates.latitude) * sin(lla_target_coordinates.latitude));

    lla_target_coordinates.latitude = radiansToDegrees(lla_target_coordinates.latitude);
    lla_target_coordinates.longitude = radiansToDegrees(lla_target_coordinates.longitude);
    lla_target_coordinates.alt = lla_start_coordinates.alt;

    double bearing_error;
    double distance_error;
    double adjusted_distance = distance;
    double adjusted_bearing = bearing;
    PolarData current_data = calculatePolarData(reference_coordinates, lla_target_coordinates);
    bool converge = false;
    double distance_error_threshold = 0.5;
    double bearing_error_threshold = 1.0 / 18.0;
    int max_loop_iterate = 100;
    int i = 0;

    while (!converge)
    {
        bearing_error = radiansToDegrees(bearing) - current_data.bearing;
        distance_error = distance - current_data.haversine;
        if (fabs(bearing_error) < bearing_error_threshold && fabs(distance_error) < distance_error_threshold)
        {
            converge = true;
        }
        else
        {
            adjusted_bearing += degreesToRadians(bearing_error);
            adjusted_distance += distance_error;

            delta = adjusted_distance / R0;
            lla_target_coordinates.latitude = asin(sin(lla_start_coordinates.latitude) * cos(delta) + cos(lla_start_coordinates.latitude) * sin(delta) * cos(adjusted_bearing));
            lla_target_coordinates.longitude = lla_start_coordinates.longitude + atan2(sin(adjusted_bearing) * sin(delta) * cos(lla_start_coordinates.latitude), cos(delta) - sin(lla_start_coordinates.latitude) * sin(lla_target_coordinates.latitude));

            lla_target_coordinates.latitude = radiansToDegrees(lla_target_coordinates.latitude);
            lla_target_coordinates.longitude = radiansToDegrees(lla_target_coordinates.longitude);

            current_data = calculatePolarData(reference_coordinates, lla_target_coordinates);
        }

        i++;
        if (i > max_loop_iterate)
        {
            converge = true;
        }
    }

    return lla_target_coordinates;
}

NavigationCalculator::LLA NavigationCalculator::findLocationFromPolarData(LLA lla_start_coordinates, double distance, double bearing)
{
    return calculateLocation(lla_start_coordinates, distance, bearing);
}

NavigationCalculator::ECEF NavigationCalculator::findLocationFromPolarData(ECEF ecef_start_coordiantes, double distance, double bearing)
{
    LLA lla_start_coodinates = convertECEF2LLA(ecef_start_coordiantes);
    return convertLLA2ECEF(calculateLocation(lla_start_coodinates, distance, bearing));
}

NavigationCalculator::LLADMS NavigationCalculator::findLocationFromPolarData(LLADMS lladms_start_coordinates, double distance, double bearing)
{
    LLA lla_start_coordinates = convertLLADMS2LLA(lladms_start_coordinates);
    return convertLLA2LLADMS(calculateLocation(lla_start_coordinates, distance, bearing));
}

NavigationCalculator::UTM NavigationCalculator::findLocationFromPolarData(UTM utm_start_coordinates, double distance, double bearing)
{
    LLA lla_start_coordinates = convertUTM2LLA(utm_start_coordinates);
    return convertLLA2UTM(calculateLocation(lla_start_coordinates, distance, bearing));
}

NavigationCalculator::MGRS NavigationCalculator::findLocationFromPolarData(MGRS mgrs_start_coordinates, double distance, double bearing)
{
    LLA lla_start_coordinates = convertMGRS2LLA(mgrs_start_coordinates);
    return convertLLA2MGRS(calculateLocation(lla_start_coordinates, distance, bearing));
}

std::string NavigationCalculator::getPolarDataString(LLA start_point, LLA target_point)
{
    PolarData data = calculatePolarData(start_point, target_point);
    return std::to_string(data.bearing) + '/' + std::to_string(data.haversine) + '/' + std::to_string(data.displacement);
}
