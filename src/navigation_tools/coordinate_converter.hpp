#ifndef COORDINATE_CONVERTER_HPP
#define COORDINATE_CONVERTER_HPP

#include <cmath>
#include <vector>
#include <iomanip>
#include <string>

#include "mgrs_lookup_functions.hpp"

class CoordinateConverter
{
public:
    // Structs for coordinate data
    struct UTM
    {
        double easting;
        double northing;
        int grid_zone;
        double alt;
    };

    struct MGRS
    {
        double easting;
        double northing;
        int grid_zone;
        double alt;
        std::string grid_letter;
        std::string false_easting;
        std::string false_northing;
    };

    struct LLA
    {
        double latitude;
        double longitude;
        double alt;
    };

    struct LLADMS
    {
        int latitude_degrees;
        int latitude_minutes;
        double latitude_seconds;
        int longitude_degrees;
        int longitude_minutes;
        double longitude_seconds;
        double alt;
    };

    struct ECEF
    {
        double x;
        double y;
        double z;
    };

    struct Vector3
    {
        double x;
        double y;
        double z;
    };

    CoordinateConverter();
    ~CoordinateConverter();

    // UTM to MGRS Conversion
    MGRS convertUTM2MGRS(UTM utm_in);

    // MGRS to UTM Conversion
    UTM convertMGRS2UTM(MGRS mgrs_in);

    // LLA to UTM Conversion
    UTM convertLLA2UTM(LLA lla_in);

    // UTM to LLA Conversion
    LLA convertUTM2LLA(UTM utm_in);

    // LLADMS to LLA
    LLA convertLLADMS2LLA(LLADMS lladms_in);

    // LLA to LLADMS
    LLADMS convertLLA2LLADMS(LLA lla_in);

    // LLA to MGRS
    MGRS convertLLA2MGRS(LLA lla_in);

    // MGRS to LLA
    LLA convertMGRS2LLA(MGRS mgrs_in);

    // ECEF to LLA
    LLA convertECEF2LLA(ECEF ecef_in);

    // LLA to ECEF
    ECEF convertLLA2ECEF(LLA lla_in);

    // ECEF TO UTM
    UTM convertECEF2UTM(ECEF ecef_in);

    // LOCAL FRAMES - ENU and NED

    double radiansToDegrees(double radians);
    double degreesToRadians(double degrees);

    // TO ENU
    Vector3 convertECEF2ENU(ECEF ecef_in, ECEF ref);
    Vector3 convertECEF2ENU(ECEF ecef_in, LLA ref);
    Vector3 convertECEF2ENU(ECEF ecef_in, LLADMS ref);
    Vector3 convertECEF2ENU(ECEF ecef_in, UTM ref);
    Vector3 convertECEF2ENU(ECEF ecef_in, MGRS ref);

    Vector3 convertLLA2ENU(LLA lla_in, ECEF ref);
    Vector3 convertLLA2ENU(LLA lla_in, LLA ref);
    Vector3 convertLLA2ENU(LLA lla_in, LLADMS ref);
    Vector3 convertLLA2ENU(LLA lla_in, UTM ref);
    Vector3 convertLLA2ENU(LLA lla_in, MGRS ref);

    Vector3 convertLLADMS2ENU(LLADMS lladms_in, ECEF ref);
    Vector3 convertLLADMS2ENU(LLADMS lladms_in, LLA ref);
    Vector3 convertLLADMS2ENU(LLADMS lladms_in, LLADMS ref);
    Vector3 convertLLADMS2ENU(LLADMS lladms_in, UTM ref);
    Vector3 convertLLADMS2ENU(LLADMS lladms_in, MGRS ref);

    Vector3 convertUTM2ENU(UTM utm_in, ECEF ref);
    Vector3 convertUTM2ENU(UTM utm_in, LLA ref);
    Vector3 convertUTM2ENU(UTM utm_in, LLADMS ref);
    Vector3 convertUTM2ENU(UTM utm_in, UTM ref);
    Vector3 convertUTM2ENU(UTM utm_in, MGRS ref);

    Vector3 convertMGRS2ENU(MGRS mgrs_in, ECEF ref);
    Vector3 convertMGRS2ENU(MGRS mgrs_in, LLA ref);
    Vector3 convertMGRS2ENU(MGRS mgrs_in, LLADMS ref);
    Vector3 convertMGRS2ENU(MGRS mgrs_in, UTM ref);
    Vector3 convertMGRS2ENU(MGRS mgrs_in, MGRS ref);

    Vector3 convertNED2ENU(Vector3 ned_in);

    // TO NED

    Vector3 convertECEF2NED(ECEF ecef_in, ECEF ref);
    Vector3 convertECEF2NED(ECEF ecef_in, LLA ref);
    Vector3 convertECEF2NED(ECEF ecef_in, LLADMS ref);
    Vector3 convertECEF2NED(ECEF ecef_in, UTM ref);
    Vector3 convertECEF2NED(ECEF ecef_in, MGRS ref);

    Vector3 convertLLA2NED(LLA lla_in, ECEF ref);
    Vector3 convertLLA2NED(LLA lla_in, LLA ref);
    Vector3 convertLLA2NED(LLA lla_in, LLADMS ref);
    Vector3 convertLLA2NED(LLA lla_in, UTM ref);
    Vector3 convertLLA2NED(LLA lla_in, MGRS ref);

    Vector3 convertLLADMS2NED(LLADMS lladms_in, ECEF ref);
    Vector3 convertLLADMS2NED(LLADMS lladms_in, LLA ref);
    Vector3 convertLLADMS2NED(LLADMS lladms_in, LLADMS ref);
    Vector3 convertLLADMS2NED(LLADMS lladms_in, UTM ref);
    Vector3 convertLLADMS2NED(LLADMS lladms_in, MGRS ref);

    Vector3 convertUTM2NED(UTM utm_in, ECEF ref);
    Vector3 convertUTM2NED(UTM utm_in, LLA ref);
    Vector3 convertUTM2NED(UTM utm_in, LLADMS ref);
    Vector3 convertUTM2NED(UTM utm_in, UTM ref);
    Vector3 convertUTM2NED(UTM utm_in, MGRS ref);

    Vector3 convertMGRS2NED(MGRS mgrs_in, ECEF ref);
    Vector3 convertMGRS2NED(MGRS mgrs_in, LLA ref);
    Vector3 convertMGRS2NED(MGRS mgrs_in, LLADMS ref);
    Vector3 convertMGRS2NED(MGRS mgrs_in, UTM ref);
    Vector3 convertMGRS2NED(MGRS mgrs_in, MGRS ref);

    Vector3 convertENU2NED(Vector3 enu_in);

    // TODO: Rotations back into global frames
    /*
    // FROM ENU
    ECEF convertENU2ECEF(Vector3 enu_in, ECEF ref);
    ECEF convertENU2ECEF(Vector3 enu_in, LLA ref);
    ECEF convertENU2ECEF(Vector3 enu_in, LLADMS ref);
    ECEF convertENU2ECEF(Vector3 enu_in, UTM ref);
    ECEF convertENU2ECEF(Vector3 enu_in, MGRS ref);

    LLA convertENU2LLA(Vector3 enu_in, ECEF ref);
    LLA convertENU2LLA(Vector3 enu_in, LLA ref);
    LLA convertENU2LLA(Vector3 enu_in, LLADMS ref);
    LLA convertENU2LLA(Vector3 enu_in, UTM ref);
    LLA convertENU2LLA(Vector3 enu_in, MGRS ref);

    LLADMS convertENU2LLADMS(Vector3 enu_in, ECEF ref);
    LLADMS convertENU2LLADMS(Vector3 enu_in, LLA ref);
    LLADMS convertENU2LLADMS(Vector3 enu_in, LLADMS ref);
    LLADMS convertENU2LLADMS(Vector3 enu_in, UTM ref);
    LLADMS convertENU2LLADMS(Vector3 enu_in, MGRS ref);

    UTM convertENU2UTM(Vector3 enu_in, ECEF ref);
    UTM convertENU2UTM(Vector3 enu_in, LLA ref);
    UTM convertENU2UTM(Vector3 enu_in, LLADMS ref);
    UTM convertENU2UTM(Vector3 enu_in, UTM ref);
    UTM convertENU2UTM(Vector3 enu_in, MGRS ref);

    MGRS convertENU2MGRS(Vector3 enu_in, ECEF ref);
    MGRS convertENU2MGRS(Vector3 enu_in, LLA ref);
    MGRS convertENU2MGRS(Vector3 enu_in, LLADMS ref);
    MGRS convertENU2MGRS(Vector3 enu_in, UTM ref);
    MGRS convertENU2MGRS(Vector3 enu_in, MGRS ref);

    // FROM NED
    ECEF convertNED2ECEF(Vector3 ned_in, ECEF ref);
    ECEF convertNED2ECEF(Vector3 ned_in, LLA ref);
    ECEF convertNED2ECEF(Vector3 ned_in, LLADMS ref);
    ECEF convertNED2ECEF(Vector3 ned_in, UTM ref);
    ECEF convertNED2ECEF(Vector3 ned_in, MGRS ref);

    LLA convertNED2LLA(Vector3 ned_in, ECEF ref);
    LLA convertNED2LLA(Vector3 ned_in, LLA ref);
    LLA convertNED2LLA(Vector3 ned_in, LLADMS ref);
    LLA convertNED2LLA(Vector3 ned_in, UTM ref);
    LLA convertNED2LLA(Vector3 ned_in, MGRS ref);

    LLADMS convertNED2LLADMS(Vector3 ned_in, ECEF ref);
    LLADMS convertNED2LLADMS(Vector3 ned_in, LLA ref);
    LLADMS convertNED2LLADMS(Vector3 ned_in, LLADMS ref);
    LLADMS convertNED2LLADMS(Vector3 ned_in, UTM ref);
    LLADMS convertNED2LLADMS(Vector3 ned_in, MGRS ref);

    UTM convertNED2UTM(Vector3 ned_in, ECEF ref);
    UTM convertNED2UTM(Vector3 ned_in, LLA ref);
    UTM convertNED2UTM(Vector3 ned_in, LLADMS ref);
    UTM convertNED2UTM(Vector3 ned_in, UTM ref);
    UTM convertNED2UTM(Vector3 ned_in, MGRS ref);

    MGRS convertNED2MGRS(Vector3 ned_in, ECEF ref);
    MGRS convertNED2MGRS(Vector3 ned_in, LLA ref);
    MGRS convertNED2MGRS(Vector3 ned_in, LLADMS ref);
    MGRS convertNED2MGRS(Vector3 ned_in, UTM ref);
    MGRS convertNED2MGRS(Vector3 ned_in, MGRS ref);
    */

private:
    // Conversion Variables
    // Earth Parameters
    double a;
    double f;
    double b;
    double e;
    double n;
    double K0;

    // IPCalculatorConstants
    double X0;
    double Y0;
    double L0;

    // Powers of e
    double e2;
    double e4;
    double e6;
    double e8;

    // Cofficients
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;

    // Other Variables
    double L;
    double z_e;
    double z_n;
    double Z_e;
    double Z_n;
    double rad2deg;
    double deg2rad;

    // Inputs and Outputs
    double lat_deg;
    double lon_deg;
    double lat_rad;
    double lon_rad;
    double easting;
    double northing;
    int grid_zone;

    // ECEF2LLA Variables
    // double rho;
    double count;
    double old;
    double error_threshold;
    bool converge;
    double N;
    double s;
    double beta;
    double den;
    double num;
    double f_inv;

    // Utility Functions
    void calculateGridZone();

    void setCoeficients1();
    void setCoeficients2();

    int sign(double x);

    bool checkGridDesignator(UTM utm_in, MGRS mgrs_in);

    Vector3 convert2ENU(ECEF ecef_point, LLA lla_reference);
};

#endif