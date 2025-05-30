#include "../include/eqnequinox.h"
#include <cmath> 

/**
 * @file eqn_equinox.cpp
 * @brief Implementación de la función para calcular la ecuación de los equinoccios.
 */

double EqnEquinox(double Mjd_TT) {
    // Obtener los ángulos de nutación (dpsi y deps)
    // NutAngles devuelve una estructura NutationAnglesResult
    NutationAnglesResult nut_angles = NutAngles(Mjd_TT);
    double dpsi = nut_angles.dpsi; // Nutación en longitud [rad]

    // Obtener la oblicuidad media de la eclíptica
    double mean_obliquity_rad = MeanObliquity(Mjd_TT); // [rad]

    // Calcular la ecuación de los equinoccios
    double EqE = dpsi * std::cos(mean_obliquity_rad);

    return EqE;
}