//
// Created by kyle on 11/26/22.
//

#ifndef mgrs_lookup_functions_hpp
#define mgrs_lookup_functions_hpp

#include <string>
#include <cmath>
// ------------------------------------------------------------------------------------------
// LOOKUP TABLES FOR MGRS
// ------------------------------------------------------------------------------------------

// UTM to MGRS
inline std::string lookupEasting(double easting_in, int grid_zone_in)
{
    int remain;
    remain = abs(grid_zone_in) % 3;

    int easting;
    easting = trunc(easting_in/100000);

    if (remain == 1)
    {
        switch (easting)
        {
            case 1:
                return "A";
                break;

            case 2:
                return "B";
                break;

            case 3:
                return "C";
                break;

            case 4:
                return "D";
                break;

            case 5:
                return "E";
                break;

            case 6:
                return "F";
                break;

            case 7:
                return "G";
                break;

            case 8:
                return "H";
                break;
        }
    }
    else if (remain == 2)
    {
        switch(easting)
        {
            case 1:
                return "J";
                break;

            case 2:
                return "K";
                break;

            case 3:
                return "L";
                break;

            case 4:
                return "M";
                break;

            case 5:
                return "N";
                break;

            case 6:
                return "P";
                break;

            case 7:
                return "Q";
                break;

            case 8:
                return "R";
                break;
        }
    }
    else if (remain == 0)
    {
        switch(easting)
        {
            case 1:
                return "S";
                break;

            case 2:
                return "T";
                break;

            case 3:
                return "U";
                break;

            case 4:
                return "V";
                break;

            case 5:
                return "W";
                break;

            case 6:
                return "X";
                break;

            case 7:
                return "Y";
                break;

            case 8:
                return "Z";
                break;
        }
    }
    return "easting not found";
}

inline std::string lookupNorthing(double northing_in, int grid_zone_in)
{
    int remain;
    remain = abs(grid_zone_in) % 2;

    int northing;
    northing = trunc(northing_in);
    northing = trunc((northing % 2000000)/100000);

    if (remain == 1)
    {
        switch(northing)
        {
            case 0:
                return "A";
                break;

            case 1:
                return "B";
                break;

            case 2:
                return "C";
                break;

            case 3:
                return "D";
                break;

            case 4:
                return "E";
                break;

            case 5:
                return "F";
                break;

            case 6:
                return "G";
                break;

            case 7:
                return "H";
                break;

            case 8:
                return "J";
                break;

            case 9:
                return "K";
                break;

            case 10:
                return "L";
                break;

            case 11:
                return "M";
                break;

            case 12:
                return "N";
                break;

            case 13:
                return "P";
                break;

            case 14:
                return "Q";
                break;

            case 15:
                return "R";
                break;

            case 16:
                return "S";
                break;

            case 17:
                return "T";
                break;

            case 18:
                return "U";
                break;

            case 19:
                return "V";
                break;
        }
    }
    else if (remain == 0)
    {
        switch(northing)
        {
            case 0:
                return "F";
                break;

            case 1:
                return "G";
                break;

            case 2:
                return "H";
                break;

            case 3:
                return "J";
                break;

            case 4:
                return "K";
                break;

            case 5:
                return "L";
                break;

            case 6:
                return "M";
                break;

            case 7:
                return "N";
                break;

            case 8:
                return "P";
                break;

            case 9:
                return "Q";
                break;

            case 10:
                return "R";
                break;

            case 11:
                return "S";
                break;

            case 12:
                return "T";
                break;

            case 13:
                return "U";
                break;

            case 14:
                return "V";
                break;

            case 15:
                return "A";
                break;

            case 16:
                return "B";
                break;

            case 17:
                return "C";
                break;

            case 18:
                return "D";
                break;

            case 19:
                return "E";
                break;
        }
    }
    return "northing not found";
}

inline std::string lookupGridZone(int latitude)
{
    if (latitude < -80)
    {
        return "NA";
    }
    else if (latitude < -72)
    {
        return "C";
    }
    else if (latitude < -64)
    {
        return "D";
    }
    else if (latitude < -56)
    {
        return "E";
    }
    else if (latitude < -48)
    {
        return "F";
    }
    else if (latitude < -40)
    {
        return "G";
    }
    else if (latitude < -32)
    {
        return "H";
    }
    else if (latitude < -24)
    {
        return "J";
    }
    else if (latitude < -16)
    {
        return "K";
    }
    else if (latitude < -8)
    {
        return "L";
    }
    else if (latitude < 0)
    {
        return "M";
    }
    else if (latitude < 8)
    {
        return "N";
    }
    else if (latitude < 16)
    {
        return "P";
    }
    else if (latitude < 24)
    {
        return "Q";
    }
    else if (latitude < 32)
    {
        return "R";
    }
    else if (latitude < 40)
    {
        return "S";
    }
    else if (latitude < 48)
    {
        return "T";
    }
    else if (latitude < 56)
    {
        return "U";
    }
    else if (latitude < 64)
    {
        return "V";
    }
    else if (latitude < 72)
    {
        return "W";
    }
    else if (latitude < 84)
    {
        return "X";
    }
    else
    {
        return "NA";
    }
}

// MGRS to UTM
inline double decodeFalseEasting(std::string false_easting)
{
    if (false_easting == "A")
    {
        return 1.0;
    }
    else if (false_easting == "B")
    {
        return 2.0;
    }
    else if (false_easting == "C")
    {
        return 3.0;
    }
    else if (false_easting == "D")
    {
        return 4.0;
    }
    else if (false_easting == "E")
    {
        return 5.0;
    }
    else if (false_easting == "F")
    {
        return 6.0;
    }
    else if (false_easting == "G")
    {
        return 7.0;
    }
    else if (false_easting == "H")
    {
        return 8.0;
    }
    else if (false_easting == "J")
    {
        return 1.0;
    }
    else if (false_easting == "K")
    {
        return 2.0;
    }
    else if (false_easting == "L")
    {
        return 3.0;
    }
    else if (false_easting == "M")
    {
        return 4.0;
    }
    else if (false_easting == "N")
    {
        return 5.0;
    }
    else if (false_easting == "P")
    {
        return 6.0;
    }
    else if (false_easting == "Q")
    {
        return 7.0;
    }
    else if (false_easting == "R")
    {
        return 8.0;
    }
    else if (false_easting == "S")
    {
        return 1.0;
    }
    else if (false_easting == "T")
    {
        return 2.0;
    }
    else if (false_easting == "U")
    {
        return 3.0;
    }
    else if (false_easting == "V")
    {
        return 4.0;
    }
    else if (false_easting == "W")
    {
        return 5.0;
    }
    else if (false_easting == "X")
    {
        return 6.0;
    }
    else if (false_easting == "Y")
    {
        return 7.0;
    }
    else if (false_easting == "Z")
    {
        return 8.0;
    }
    return 0;
}

inline double decodeFalseNorthing(std::string false_northing, int grid_zone)
{
    int remain = grid_zone % 2;
    if (remain==1)
    {
        if (false_northing=="A")
        {
            return 0.0;
        }
        else if(false_northing=="B")
        {
            return 1.0;
        }
        else if(false_northing=="C")
        {
            return 2.0;
        }
        else if(false_northing=="D")
        {
            return 3.0;
        }
        else if(false_northing=="E")
        {
            return 4.0;
        }
        else if(false_northing=="F")
        {
            return 5.0;
        }
        else if(false_northing=="G")
        {
            return 6.0;
        }
        else if(false_northing=="H")
        {
            return 7.0;
        }
        else if(false_northing=="J")
        {
            return 8.0;
        }
        else if(false_northing=="K")
        {
            return 9.0;
        }
        else if(false_northing=="L")
        {
            return 10.0;
        }
        else if(false_northing=="M")
        {
            return 11.0;
        }
        else if(false_northing=="N")
        {
            return 12.0;
        }
        else if(false_northing=="P")
        {
            return 13.0;
        }
        else if(false_northing=="Q")
        {
            return 14.0;
        }
        else if(false_northing=="R")
        {
            return 15.0;
        }
        else if(false_northing=="S")
        {
            return 16.0;
        }
        else if(false_northing=="T")
        {
            return 17.0;
        }
        else if(false_northing=="U")
        {
            return 18.0;
        }
        else if(false_northing=="V")
        {
            return 19.0;
        }
    }
    else
    {
        if (false_northing=="F")
        {
            return 0.0;
        }
        else if(false_northing=="G")
        {
            return 1.0;
        }
        else if(false_northing=="H")
        {
            return 2.0;
        }
        else if(false_northing=="J")
        {
            return 3.0;
        }
        else if(false_northing=="K")
        {
            return 4.0;
        }
        else if(false_northing=="L")
        {
            return 5.0;
        }
        else if(false_northing=="M")
        {
            return 6.0;
        }
        else if(false_northing=="N")
        {
            return 7.0;
        }
        else if(false_northing=="P")
        {
            return 8.0;
        }
        else if(false_northing=="Q")
        {
            return 9.0;
        }
        else if(false_northing=="R")
        {
            return 10.0;
        }
        else if(false_northing=="S")
        {
            return 11.0;
        }
        else if(false_northing=="T")
        {
            return 12.0;
        }
        else if(false_northing=="U")
        {
            return 13.0;
        }
        else if(false_northing=="V")
        {
            return 14.0;
        }
        else if(false_northing=="A")
        {
            return 15.0;
        }
        else if(false_northing=="B")
        {
            return 16.0;
        }
        else if(false_northing=="C")
        {
            return 17.0;
        }
        else if(false_northing=="D")
        {
            return 18.0;
        }
        else if(false_northing=="E")
        {
            return 19.0;
        }
    }
    return 0;
}

#endif //POLARCOORDINATES_MGRS_LOOKUP_FUNCTIONS_H
