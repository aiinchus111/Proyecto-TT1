
#include "../include/r_z.h"
#include <cmath> 

/**
 * @file r_z.cpp
 * @brief Implementación de la función para generar una matriz de rotación alrededor del eje Z.
 */

Matrix R_z(double angle) {
    double C = std::cos(angle);
    double S = std::sin(angle);

    Matrix rotmat(3, 3); 

    // Fila 1
    rotmat(1, 1) =  C;   rotmat(1, 2) = S;   rotmat(1, 3) = 0.0;
    // Fila 2
    rotmat(2, 1) = -S;   rotmat(2, 2) = C;   rotmat(2, 3) = 0.0;
    // Fila 3
    rotmat(3, 1) = 0.0;  rotmat(3, 2) = 0.0; rotmat(3, 3) = 1.0;

    return rotmat;
}