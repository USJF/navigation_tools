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
    ~NavigationCalculator();

    struct PolarData
    {
        double displacement;
        double haversine;
        double bearing;
    };

    PolarData calculatePolarData(LLA start_point, LLA target_point);

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

    double radiansToDegrees(double radians);
    double degreesToRadians(double degrees);
    double calculateBearing(double delta_east, double delta_north);
};

#endif