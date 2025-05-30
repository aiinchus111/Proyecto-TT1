#include "../include/gha_matrix.h"


/**
 * @file gha_matrix.cpp
 * @brief Implementación de la función para calcular la matriz del Ángulo Horario de Greenwich.
 */

Matrix GHAMatrix(double Mjd_UT1) {
    // Calcular el Tiempo Sidéreo Aparente de Greenwich (GAST)
    double gast_rad = gast(Mjd_UT1); // [rad]

    // Construir la matriz de rotación alrededor del eje Z usando GAST
    // GHAmat = R_z(gast_rad)
    Matrix GHAmat = R_z(gast_rad);
    
    return GHAmat;
}