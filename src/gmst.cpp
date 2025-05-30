#include "../include/gmst.h"
#include <cmath> // Para std::floor

/**
 * @file gmst.cpp
 * @brief Implementación de la función para calcular el Tiempo Sidéreo Medio de Greenwich (GMST).
 */

double gmst(double Mjd_UT1) {
    const double Secs_per_day = 86400.0;

    // Calcular Mjd a las 0h UT1 del día actual
    double Mjd_0 = std::floor(Mjd_UT1);
    
    // UT1 del día actual en segundos desde las 0h UT1
    double UT1_seconds = Secs_per_day * (Mjd_UT1 - Mjd_0);
    
    // Siglos Julianos desde J2000.0 hasta las 0h UT1 del día actual
    double T_0 = (Mjd_0 - Const::MJD_J2000) / 36525.0;
    
    // Siglos Julianos desde J2000.0 hasta el instante Mjd_UT1 exacto
    double T = (Mjd_UT1 - Const::MJD_J2000) / 36525.0;

    // Fórmula para GMST en segundos sidéreos (IAU 1982/1994)
    double gmst_seconds = 24110.54841 +
                          8640184.812866 * T_0 +
                          1.002737909350795 * UT1_seconds +
                          (0.093104 - 6.2e-6 * T) * T * T;
    
    // Normalizar a [0, Secs_per_day) y luego a [0, 1) fracción de día/revolución
    double gmst_frac_of_day = Frac(gmst_seconds / Secs_per_day);
    
    // Convertir a radianes [0, 2*pi)
    double gmst_radians = Const::pi2 * gmst_frac_of_day; // Const::pi2 es 2*pi

    return gmst_radians;
}