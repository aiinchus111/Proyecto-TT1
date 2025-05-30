#include "../include/r_x.h"
#include <cmath> 

/**
 * @file r_x.cpp
 * @brief Implementación de la función para generar una matriz de rotación alrededor del eje X.
 */

Matrix R_x(double angle) {
    double C = std::cos(angle);
    double S = std::sin(angle);

    Matrix rotmat(3, 3);

    // Fila 1
    rotmat(1, 1) = 1.0; rotmat(1, 2) = 0.0; rotmat(1, 3) = 0.0;
    // Fila 2
    rotmat(2, 1) = 0.0; rotmat(2, 2) = C;   rotmat(2, 3) = S;
    // Fila 3
    rotmat(3, 1) = 0.0; rotmat(3, 2) = -S;  rotmat(3, 3) = C;

    return rotmat;
}
