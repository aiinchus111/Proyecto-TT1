#include "../include/mjday_tdb.h"
#include <cmath> 

/**
 * @file ephemerides_utils.cpp
 * @brief Implementación de utilidades para cálculos de efemérides.
 */

double Mjday_TDB(double Mjd_TT) {
    const double SECONDS_PER_DAY = 86400.0;

    // Calcular siglos Julianos de TT desde J2000.0
    double T_TT = (Mjd_TT - Const::MJD_J2000) / 36525.0;

    
    double delta_T_seconds = 
        0.001658 * std::sin(628.3076  * T_TT + 6.2401) +
        0.000022 * std::sin(575.3385  * T_TT + 4.2970) +
        0.000014 * std::sin(1256.6152 * T_TT + 6.1969) +
        0.000005 * std::sin(606.9777  * T_TT + 4.0212) +
        0.000005 * std::sin(52.9691   * T_TT + 0.4444) +
        0.000002 * std::sin(21.3299   * T_TT + 5.5431) +
        0.000010 * std::sin(628.3076  * T_TT + 4.2490); 
    double Mjd_TDB = Mjd_TT + delta_T_seconds / SECONDS_PER_DAY;

    return Mjd_TDB;
}