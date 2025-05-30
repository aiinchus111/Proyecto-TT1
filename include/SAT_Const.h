#ifndef _CONST_
#define _CONST_

#include <cmath>
#include <iostream>

#define M_PI 3.14159265358979323846


namespace Const {
    constexpr double pi = M_PI;
    
    constexpr double pi2 = 2 * M_PI;           
    constexpr double Rad = M_PI / 180.0;     
    constexpr double Deg = 180.0 / M_PI;     
    constexpr double Arcs = 3600.0 * 180.0 / M_PI; 

    constexpr double MJD_J2000 = 51544.5;      
    constexpr double T_B1950 = -0.500002108;   
    constexpr double c_light = 299792458.0;    
    constexpr double AU = 149597870700.0;       

    constexpr double R_Earth = 6378.1363e3;    
    constexpr double f_Earth = 1.0 / 298.257223563; 
    constexpr double R_Sun = 696000e3;        
    constexpr double R_Moon = 1738e3;         

    constexpr double omega_Earth = 15.04106717866910 / 3600.0 * Rad; // [rad/s]; WGS-84

    constexpr double GM_Earth = 398600.435436e9;       // [m^3/s^2]; DE430
    constexpr double GM_Sun = 132712440041.939400e9;    // [m^3/s^2]; DE430
    constexpr double GM_Moon = GM_Earth / 81.30056907419062; // [m^3/s^2]; DE430
    constexpr double GM_Mercury = 22031.780000e9;        // [m^3/s^2]; DE430
    constexpr double GM_Venus = 324858.592000e9;         // [m^3/s^2]; DE430
    constexpr double GM_Mars = 42828.375214e9;         // [m^3/s^2]; DE430
    constexpr double GM_Jupiter = 126712764.800000e9;     // [m^3/s^2]; DE430
    constexpr double GM_Saturn = 37940585.200000e9;      // [m^3/s^2]; DE430
    constexpr double GM_Uranus = 5794548.600000e9;       // [m^3/s^2]; DE430
    constexpr double GM_Neptune = 6836527.100580e9;      // [m^3/s^2]; DE430
    constexpr double GM_Pluto = 977.0000000000009e9;      // [m^3/s^2]; DE430

    constexpr double P_Sol = 1367.0 / c_light;    // [N/m^2] (~1367 W/m^2); IERS 96
}

#endif