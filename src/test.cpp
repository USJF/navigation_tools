#include <iostream>
#include <iomanip>

#include "navigation_calculator.hpp"

int main(int argc, char **argv)
{
    NavigationCalculator calculator;

    NavigationCalculator::LLA start_point;
    start_point.latitude = 34.862976;
    start_point.longitude = -86.778458;
    start_point.alt = 135;

    NavigationCalculator::LLA target_point;
    target_point.latitude = 55.042636;
    target_point.longitude = -106.983825;
    target_point.alt = 135;//202;

    NavigationCalculator::PolarData data = calculator.calculatePolarData(start_point, target_point);

    std::cout<< data.bearing << "\t" << data.haversine << std::endl;
    
    NavigationCalculator::LLA lla_output = calculator.findLocationFromPolarData(start_point, data.haversine, calculator.degreesToRadians(data.bearing));
    NavigationCalculator::PolarData data2 = calculator.calculatePolarData(start_point, lla_output);
    std::cout << lla_output.latitude << '\t' << lla_output.longitude << std::endl;
    std::cout << data2.bearing << '\t' << data2.haversine <<std::endl;

    bool converge = false;
    double distance_threshold = 1.0;
    double bearing_threshold = 1.0/18.0;
    double distance_error; 
    double bearing_error;
    int i = 0;

    NavigationCalculator::LLA point;
    point.latitude = 0;
    point.longitude = 0;
    point.alt = 0;
    calculator.setLocationLLA(point);
    std::string data_string = calculator.getLocationData();
    std::cout<<data_string<<std::endl;
    /*NavigationCalculator::PolarData converge_data_buffer = data;
    while (!converge)
    {
        //std::cout << "Pass " << i << "--------------" << std::endl;
        bearing_error = data.bearing - data2.bearing;
        distance_error = data.haversine - data2.haversine;
        //std::cout << bearing_error << '\t' << distance_error << std::endl;
        if (fabs(bearing_error) < bearing_threshold && fabs(distance_error) < distance_threshold)
        {
            converge = true;
        }
        else
        {
            converge_data_buffer.bearing = converge_data_buffer.bearing + bearing_error;
            converge_data_buffer.haversine = converge_data_buffer.haversine + distance_error;
            lla_output = calculator.findLocationFromPolarData(start_point, converge_data_buffer.haversine, calculator.degreesToRadians(converge_data_buffer.bearing));
            data2 = calculator.calculatePolarData(start_point, lla_output);
        }
        std::cout << converge_data_buffer.bearing << '\t' << converge_data_buffer.haversine << '\t' << converge << std::endl;
        i++;
        if (i > 5)
        {
            converge = true;
        }
    }
    std::cout<<i<<std::endl;
    std::cout<<lla_output.latitude << '\t' << lla_output.longitude <<std::endl;
    std::cout<<data2.bearing << '\t' << data2.haversine<<std::endl;*/
    /*
    NavigationCalculator::LLADMS lladms_output;
    lladms_output.latitude_degrees = 35;
    lladms_output.latitude_minutes = 02;
    lladms_output.latitude_seconds = 33;
    lladms_output.longitude_degrees = -86;
    lladms_output.longitude_minutes = -59;
    lladms_output.longitude_seconds = -04;
    NavigationCalculator::LLA lla_internet_answer = calculator.convertLLADMS2LLA(lladms_output);
    std::cout << lla_internet_answer.latitude << '\t' <<lla_internet_answer.longitude << std::endl;*/
}