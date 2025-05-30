#include "../include/r_y.h"
#include <cmath> 

/**
 * @file r_y.cpp
 * @brief Implementación de la función para generar una matriz de rotación alrededor del eje Y.
 */

Matrix R_y(double angle) {
    double C = std::cos(angle);
    double S = std::sin(angle);

    Matrix rotmat(3, 3); 

    // Fila 1
    rotmat(1, 1) =  C;   rotmat(1, 2) = 0.0; rotmat(1, 3) = -S;
    // Fila 2
    rotmat(2, 1) = 0.0;  rotmat(2, 2) = 1.0; rotmat(2, 3) = 0.0;
    // Fila 3
    rotmat(3, 1) =  S;   rotmat(3, 2) = 0.0; rotmat(3, 3) =  C;

    return rotmat;
}