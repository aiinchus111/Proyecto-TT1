#include "../include/polematrix.h"


/**
 * @file polematrix.cpp
 * @brief Implementación de la función para calcular la matriz de movimiento polar.
 */

Matrix PoleMatrix(double xp, double yp) {
    // Calcular las matrices de rotación elemental con los ángulos negativos
    Matrix Rx_neg_yp = R_x(-yp);
    Matrix Ry_neg_xp = R_y(-xp);


    Matrix PoleMat = Ry_neg_xp * Rx_neg_yp;
    
    return PoleMat;
}