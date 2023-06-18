#include "coordinate_converter.hpp"

CoordinateConverter::CoordinateConverter()
{
    a = 6378137.0;
    f = 298.257223563;
    b = a * (1.0 - 1.0 / f);
    K0 = 0.9996;
    X0 = 500000;

    e = sqrt((a * a - b * b) / (a * a));

    e2 = pow(e, 2);
    e4 = pow(e, 4);
    e6 = pow(e, 6);
    e8 = pow(e, 8);

    n = a * K0;

    error_threshold = 1e-16;
    rad2deg = 180.0 / M_PI;
    deg2rad = M_PI / 180.0;

    // ECEF vars
    f_inv = 1.0/f;
}

CoordinateConverter::~CoordinateConverter()
{
    // nothing
}

void CoordinateConverter::calculateGridZone()
{
    int bound = -174;
    int i = 1;
    bool key = true;

    // Calculate the "Nominal" Grid Zone
    while (key)
    {
        if (lon_deg < bound)
        {
            grid_zone = i;
            key = false;
        }
        if (bound > 180)
        {
            key = false;
        }
        i++;
        bound += 6;
    }

    // Administrative Exceptions to Grid Zone
    if ((grid_zone == 31 && lon_deg >= 3.0) && (lat_deg < 64.0 && lat_deg >= 56.0))
    {
        grid_zone = 32;
        L0 = 9;
    }
    else if (grid_zone == 32 && lat_deg >= 72.0)
    {
        if (lon_deg < 9.0)
        {
            grid_zone = 31;
        }
        else
        {
            grid_zone = 33;
        }
    }
    else if (grid_zone == 34 && lat_deg >= 72.0)
    {
        if (lon_deg < 21.0)
        {
            grid_zone = 33;
        }
        else
        {
            grid_zone = 35;
        }
    }
    else if (grid_zone == 36 && lat_deg >= 72.0)
    {
        if (lon_deg < 33.0)
        {
            grid_zone = 35;
        }
        else
        {
            grid_zone = 37;
        }
    }

    // Calculate the Central Meridian for the Grid Zone
    L0 = (grid_zone * 6.0) - 183.0;
}

void CoordinateConverter::setCoeficients1()
{
    c1 = (-175.0 / 16384.0) * e8 + (-5.0 / 256.0) * e6 + (-3.0 / 64.0) * e4 + (-1.0 / 4.0) * e2 + 1;
    c2 = (-901.0 / 184320.0) * e8 + (-9.0 / 1024.0) * e6 + (-1.0 / 96.0) * e4 + (1.0 / 8.0) * e2;
    c3 = (-311.0 / 737280.0) * e8 + (17.0 / 5120.0) * e6 + (13.0 / 768.0) * e4;
    c4 = (899.0 / 430080.0) * e8 + (61.0 / 15360.0) * e6;
    c5 = (49561.0 / 41287680.0) * e8;
}

void CoordinateConverter::setCoeficients2()
{
    c1 = (-175.0 / 16384.0) * e8 + (-5.0 / 256.0) * e6 + (-3.0 / 64.0) * e4 + (-1.0 / 4.0) * e2 + 1;
    c2 = (1.0 / 61440.0) * e8 + (7.0 / 2048.0) * e6 + (1.0 / 48.0) * e4 + (1.0 / 8.0) * e2;
    c3 = (559.0 / 368640.0) * e8 + (3.0 / 1280.0) * e6 + (1.0 / 768.0) * e4;
    c4 = (283.0 / 430080.0) * e8 + (17.0 / 30720.0) * e6;
    c5 = (4397.0 / 41287680.0) * e8;
}

int CoordinateConverter::sign(double x)
{
    if (x >= 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

bool CoordinateConverter::checkGridDesignator(CoordinateConverter::UTM utm_in, CoordinateConverter::MGRS mgrs_in)
{
    MGRS check_coords = convertUTM2MGRS(utm_in);

    if (check_coords.grid_letter != mgrs_in.grid_letter)
    {
        return false;
    }

    if (check_coords.false_northing != mgrs_in.false_northing)
    {
        return false;
    }

    if (check_coords.false_easting != mgrs_in.false_easting)
    {
        return false;
    }

    return true;
}

CoordinateConverter::MGRS CoordinateConverter::convertUTM2MGRS(CoordinateConverter::UTM utm_in)
{
    int false_easting;
    int false_northing;

    // Get LLA Coordinates for the grid zone letter
    LLA lla_coords = convertUTM2LLA(utm_in);

    false_easting = trunc(utm_in.easting / 100000);
    false_northing = trunc(utm_in.northing / 100000);

    MGRS mgrs_coordinates;
    mgrs_coordinates.easting = utm_in.easting - false_easting * 100000;
    mgrs_coordinates.northing = utm_in.northing - false_northing * 100000;
    mgrs_coordinates.false_easting = lookupEasting(utm_in.easting, utm_in.grid_zone);
    mgrs_coordinates.false_northing = lookupNorthing(utm_in.northing, utm_in.grid_zone);
    mgrs_coordinates.grid_zone = utm_in.grid_zone;
    mgrs_coordinates.grid_letter = lookupGridZone(lla_coords.latitude);
    mgrs_coordinates.alt = utm_in.alt;

    return mgrs_coordinates;
}

CoordinateConverter::UTM CoordinateConverter::convertMGRS2UTM(CoordinateConverter::MGRS mgrs_in)
{
    UTM utm_coords;
    utm_coords.grid_zone = mgrs_in.grid_zone;
    utm_coords.easting = mgrs_in.easting + decodeFalseEasting(mgrs_in.false_easting) * 100000.0;
    utm_coords.northing = mgrs_in.northing + decodeFalseNorthing(mgrs_in.false_northing, mgrs_in.grid_zone) * 100000.0;
    utm_coords.alt = mgrs_in.alt;

    while (!checkGridDesignator(utm_coords, mgrs_in) && utm_coords.northing < 10000000)
    {
        utm_coords.northing += 2000000;
    }

    return utm_coords;
}

CoordinateConverter::UTM CoordinateConverter::convertLLA2UTM(CoordinateConverter::LLA lla_in)
{
    lat_deg = lla_in.latitude;
    lon_deg = lla_in.longitude;

    setCoeficients1();

    if (lat_deg < 0)
    {
        Y0 = 10000000.0;
    }
    else
    {
        Y0 = 0;
    }

    calculateGridZone();

    lon_rad = lon_deg * deg2rad;
    lat_rad = lat_deg * deg2rad;
    L0 = L0 * deg2rad;

    L = log(tan(M_PI / 4 + lat_rad / 2) * (pow(((1 - e * sin(lat_rad)) / (1 + e * sin(lat_rad))), (e / 2))));
    z_e = log(tan(M_PI / 4 + asin(sin(lon_rad - L0) / cosh(L)) / 2));
    z_n = atan(sinh(L) / cos(lon_rad - L0));

    Z_e = n * c1 * z_e + n * (c2 * cos(2 * z_n) * sinh(2 * z_e) + c3 * cos(4 * z_n) * sinh(4 * z_e) + c4 * cos(6 * z_n) * sinh(6 * z_e) + c5 * cos(8 * z_n) * sinh(8 * z_e));
    Z_n = n * c1 * z_n + n * (c2 * sin(2 * z_n) * cosh(2 * z_e) + c3 * sin(4 * z_n) * cosh(4 * z_e) + c4 * sin(6 * z_n) * cosh(6 * z_e) + c5 * sin(8 * z_n) * cosh(8 * z_e));

    easting = Z_e + X0;
    northing = Z_n + Y0;

    UTM utm_coordinates;
    utm_coordinates.easting = easting;
    utm_coordinates.northing = northing;
    if (lat_deg >= 0)
    {
        utm_coordinates.grid_zone = grid_zone;
    }
    else
    {
        utm_coordinates.grid_zone = -1 * grid_zone;
    }
    utm_coordinates.alt = lla_in.alt;

    return utm_coordinates;
}

CoordinateConverter::LLA CoordinateConverter::convertUTM2LLA(CoordinateConverter::UTM utm_in)
{
    easting = utm_in.easting;
    northing = utm_in.northing;
    grid_zone = utm_in.grid_zone;

    // Hemisphere Based Northing Offset
    if (grid_zone >= 0.0)
    {
        Y0 = 0;
    }
    else
    {
        Y0 = 10000000.0;
    }

    // Central Meridian
    L0 = (6.0 * abs(grid_zone) - 183.0) * deg2rad;

    // Parameters for converging equation
    int max_nubmer_iterations;
    max_nubmer_iterations = 10000;
    double convergence_threshold;
    convergence_threshold = pow(10.0, -100.0);

    // Set Coefficients
    setCoeficients2();

    // Calculate Footpring LLA
    z_n = (northing - Y0) / n / c1;
    z_e = (easting - X0) / n / c1;

    Z_n = z_n - c2 * sin(2.0 * z_n) * cosh(2.0 * z_e) - c3 * sin(4.0 * z_n) * cosh(4.0 * z_e) - c4 * sin(6.0 * z_n) * cosh(6.0 * z_e) - c5 * sin(8.0 * z_n) * cosh(8.0 * z_e);
    Z_e = z_e - c2 * cos(2.0 * z_n) * sinh(2.0 * z_e) - c3 * cos(4.0 * z_n) * sinh(4.0 * z_e) - c4 * cos(6.0 * z_n) * sinh(6.0 * z_e) - c5 * cos(8.0 * z_n) * sinh(8.0 * z_e);

    double p;

    p = asin(sin(Z_n) / cosh(Z_e));
    L = log(tan(M_PI / 4.0 + p / 2.0));

    lat_rad = 2.0 * atan(exp(L)) - M_PI / 2.0;
    lon_rad = L0 + atan(sinh(Z_e) / cos(Z_n));

    // Converge Latitude
    double lat0;
    lat0 = 0;

    int i = 0;
    double es;

    while (sqrt(pow(lat_rad - lat0, 2)) > convergence_threshold && i < max_nubmer_iterations)
    {
        lat0 = lat_rad;
        es = e * sin(lat0);
        lat_rad = 2 * atan(pow(((1.0 + es) / (1.0 - es)), (e / 2.0)) * exp(L)) - M_PI / 2.0;
        i++;
    }

    lat_deg = lat_rad * rad2deg;
    lon_deg = lon_rad * rad2deg;

    LLA lla_coords;
    lla_coords.latitude = lat_deg;
    lla_coords.longitude = lon_deg;
    lla_coords.alt = utm_in.alt;

    return lla_coords;
}

CoordinateConverter::LLA CoordinateConverter::convertLLADMS2LLA(CoordinateConverter::LLADMS lladms_in)
{
    LLA lla_coords;
    lla_coords.latitude = lladms_in.latitude_degrees + lladms_in.latitude_minutes / 60.0 + lladms_in.latitude_seconds / (60.0 * 60.0);
    lla_coords.longitude = lladms_in.longitude_degrees + lladms_in.longitude_minutes / 60.0 + lladms_in.longitude_seconds / (60.0 * 60.0);
    lla_coords.alt = lladms_in.alt;
    return lla_coords;
}

CoordinateConverter::LLADMS CoordinateConverter::convertLLA2LLADMS(CoordinateConverter::LLA lla_in)
{
    LLADMS lladms_coords;

    lladms_coords.latitude_degrees = trunc(lla_in.latitude);
    double remain = lla_in.latitude - lladms_coords.latitude_degrees;
    lladms_coords.latitude_minutes = trunc(remain * 60.0);
    remain = remain - lladms_coords.latitude_minutes / 60.0;
    lladms_coords.latitude_seconds = remain * (60.0 * 60.0);

    lladms_coords.longitude_degrees = trunc(lla_in.longitude);
    remain = lla_in.longitude - lladms_coords.longitude_degrees;
    lladms_coords.longitude_minutes = trunc(remain * 60.0);
    remain = remain - lladms_coords.longitude_minutes / 60.0;
    lladms_coords.longitude_seconds = remain * (60.0 * 60.0);

    lladms_coords.alt = lla_in.alt;

    return lladms_coords;
}

CoordinateConverter::MGRS CoordinateConverter::convertLLA2MGRS(CoordinateConverter::LLA lla_in)
{
    return convertUTM2MGRS(convertLLA2UTM(lla_in));
}

CoordinateConverter::LLA CoordinateConverter::convertMGRS2LLA(CoordinateConverter::MGRS mgrs_in)
{
    return convertUTM2LLA(convertMGRS2UTM(mgrs_in));
}

CoordinateConverter::LLA CoordinateConverter::convertECEF2LLA(CoordinateConverter::ECEF ecef_in)
{
    LLA lla_coords;
    lla_coords.longitude = atan2(ecef_in.y, ecef_in.x) * rad2deg;

    s = sqrt(pow(ecef_in.x, 2) + pow(ecef_in.y, 2));

    beta = atan2(ecef_in.z, (1-f_inv)*s);
    num = ecef_in.z + (e2 * (1-f_inv)/(1-e2)) * a * pow(sin(beta), 3);
    den = s - e2 * a * pow(cos(beta), 3);
    lat_rad = atan2(num, den);

    converge = false;
    count = 0;
    while (!converge)
    {
        old = lat_rad;
        beta = atan2((1-f_inv)*sin(lat_rad), cos(lat_rad));
        num = ecef_in.z + (e2 * (1-f_inv)/(1-e2)) * a * pow(sin(beta), 3);
        den = s - e2 * a * pow(cos(beta), 3);
        lat_rad = atan2(num, den);

        if (abs(lat_rad - old) * rad2deg < error_threshold)
        {
            converge = true;
        }   
        else if(count > 100)
        {
            converge = true;
        }
        count++;
    }

    lla_coords.latitude = lat_rad*rad2deg;
    N = a/sqrt(1-e2*pow(sin(lat_rad), 2));
    lla_coords.alt = s * cos(lat_rad) + (ecef_in.z + e2 * N * sin(lat_rad)) * sin(lat_rad) - N;

    return lla_coords;
}

CoordinateConverter::ECEF CoordinateConverter::convertLLA2ECEF(CoordinateConverter::LLA lla_in)
{
    ECEF ecef_coords;

    N = a / sqrt(1 - e2 * pow(sin(lla_in.latitude * deg2rad), 2));

    ecef_coords.x = (N + lla_in.alt) * cos(lla_in.latitude * deg2rad) * cos(lla_in.longitude * deg2rad);
    ecef_coords.y = (N + lla_in.alt) * cos(lla_in.latitude * deg2rad) * sin(lla_in.longitude * deg2rad);
    ecef_coords.z = ((1-e2)*N + lla_in.alt) * sin(lla_in.latitude * deg2rad);

    return ecef_coords;
}

CoordinateConverter::UTM CoordinateConverter::convertECEF2UTM(CoordinateConverter::ECEF ecef_in)
{
    LLA lla_coords;
    lla_coords = convertECEF2LLA(ecef_in);

    return convertLLA2UTM(lla_coords);
}

