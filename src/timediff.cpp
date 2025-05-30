#include "../include/timediff.h"

/**
 * @file timediff.cpp
 * @brief Implementación de la función para calcular diferencias entre escalas de tiempo.
 */

TimeDifferences timediff(double UT1_UTC, double TAI_UTC) {
    // Constantes internas definidas en la función original de MATLAB
    const double TT_TAI  = +32.184; // TT - TAI [s]
    const double GPS_TAI = -19.0;   // GPS Time - TAI [s]

    TimeDifferences results;

    
    double UTC_TAI = -TAI_UTC;         // UTC - TAI

    // Cálculos de las diferencias a devolver
    results.UT1_TAI = UT1_UTC - TAI_UTC;
    results.UTC_GPS = UTC_TAI - GPS_TAI;
    results.UT1_GPS = results.UT1_TAI - GPS_TAI; 
    results.TT_UTC  = TT_TAI - UTC_TAI;
    results.GPS_UTC = GPS_TAI - UTC_TAI;

    return results;
}