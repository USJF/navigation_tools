#ifndef NAVIGATION_CALCULATOR_HPP
#define NAVIGATION_CALCULATOR_HPP

#include "coordinate_converter.hpp"
#include <vector>
#include <cmath>
#include <string>

class NavigationCalculator : public CoordinateConverter
{
    public:
    NavigationCalculator();
    NavigationCalculator(LLA lla_initial_position);
    NavigationCalculator(ECEF ecef_initial_position);
    NavigationCalculator(UTM utm_initial_position);
    NavigationCalculator(LLADMS lladms_initial_position);
    NavigationCalculator(MGRS mgrs_initial_position);

    ~NavigationCalculator();

    struct PolarData
    {
        double displacement;
        double haversine;
        double bearing;
    };

    struct Vector3
    {
        double x;
        double y;
        double z;
    };

    PolarData calculatePolarData(LLA start_point, LLA target_point);
    double radiansToDegrees(double radians);
    double degreesToRadians(double degrees);

    std::string output_string;
    std::string getLocationData();

    double calculateRn(double latitude);
    double calculateRe(double latitude);
    double haversine(double latitude1, double latitude2, double longitude1, double longitude2);

    void setLocationECEF(double x, double y, double z);
    void setLocationECEF(ECEF ecef_coordinates);

    void setLocationUTM(double easting, double northing, double alt, int grid_zone);
    void setLocationUTM(UTM utm_coordinates);

    void setLocationMGRS(double easting, double northing, double alt, int grid_number, std::string grid_letter, std::string false_easting, std::string false_northing);
    void setLocationMGRS(double easting, double northing, double alt, int grid_number, std::string grid_letters);
    void setLocationMGRS(MGRS mgrs_coordinates);

    void setLocationLLA(double latitude, double longitude, double alt);
    void setLocationLLA(LLA lla_coordinates);
    
    void setLocationLLADMS(int latitude_degrees, int latitude_minutes, double latitude_seconds, int longitude_degrees, int longitude_minutes, double longitude_seconds, double alt);
    void setLocationLLADMS(LLADMS lladms_coordinates);

    Vector3 getBodyRPV(LLA lla_target_coordinates, double heading);
    Vector3 getBodyRPV(LLA lla_target_coordinates, double heading, LLA lla_start_coordinates);
    Vector3 getBodyRPV(ECEF ecef_target_coordinates, double heading);
    Vector3 getBodyRPV(ECEF ecef_target_coordinates, double heading, ECEF ecef_start_coordinates);
    Vector3 getBodyRPV(LLADMS lladms_target_coordinates, double heading);
    Vector3 getBodyRPV(LLADMS lladms_target_coordinates, double heading, LLADMS lladms_start_coordinates);
    Vector3 getBodyRPV(UTM utm_target_coordinates, double heading);
    Vector3 getBodyRPV(UTM utm_target_coordinates, double heading, UTM utm_start_coordinates);
    Vector3 getBodyRPV(MGRS mgrs_target_coordinates, double heading);
    Vector3 getBodyRPV(MGRS mgrs_target_coordinates, double heading, MGRS mgrs_start_coordinates);

    LLA calculateLocation(LLA lla_start_coordinates, double distance, double bearing);

    LLA findLocationFromPolarData(LLA lla_start_coordinates, double distance, double bearing);
    ECEF findLocationFromPolarData(ECEF ecef_start_coordinates, double distance, double bearing);
    LLADMS findLocationFromPolarData(LLADMS lladms_start_coordinates, double distance, double bearing);
    UTM findLocationFromPolarData(UTM utm_start_coordinates, double distance, double bearing);
    MGRS findLocationFromPolarData(MGRS mgrs_start_coordinates, double distance, double bearing);

    private:
    double Rp;
    double R0;
    double e;
    double e2;

    double delta_x;
    double delta_y;
    double delta_z;

    ECEF location_ecef;
    UTM location_utm;
    MGRS location_mgrs;
    LLA location_lla;
    LLADMS location_lladms;

    // Polar Calculation Variables 
    PolarData output_data;
    double local_rotation_coefs[3][3];
    ECEF ecef_start;
    ECEF ecef_target;
    double latitude;
    double longitude;
    double magnitude;
    double enu[3];

    double calculateBearing(double delta_east, double delta_north);

    Vector3 calculateLocalRPV(LLA start_point, LLA target_point);
    Vector3 calculateLocalRPV(ECEF start_point, ECEF target_point);
    Vector3 calculateLocalRPV(LLADMS start_point, LLADMS target_point);
    Vector3 calculateLocalRPV(UTM start_point, UTM target_point);
    Vector3 calculateLocalRPV(MGRS start_point, MGRS target_point);
};

#endif