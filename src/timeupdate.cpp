#include "../include/timeupdate.h"

/**
 * @file time_update.cpp
 * @brief Implementación de las funciones para la actualización temporal de la matriz de covarianza.
 */

Matrix TimeUpdate( Matrix& P_current,  Matrix& Phi,  Matrix& Qdt) {
    // P_new = Phi * P_current * Phi_transpose + Qdt
    
    Matrix Phi_transpose = Phi.transpose();
    Matrix P_new = Phi * P_current * Phi_transpose + Qdt;
    
    return P_new;
}

Matrix TimeUpdate( Matrix& P_current,  Matrix& Phi) {

    
    Matrix Phi_transpose = Phi.transpose();
    Matrix P_new = Phi * P_current * Phi_transpose;
    
    return P_new;
}