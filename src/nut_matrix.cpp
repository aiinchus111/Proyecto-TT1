#include "../include/nut_matrix.h"


/**
 * @file nut_matrix.cpp
 * @brief Implementación de la función para calcular la matriz de nutación.
 */

Matrix NutMatrix(double Mjd_TT) {
    // Oblicuidad media de la eclíptica
    double eps_mean = MeanObliquity(Mjd_TT); // [rad]

    // Ángulos de nutación en longitud y oblicuidad
    NutationAnglesResult nut_angles = NutAngles(Mjd_TT);
    double dpsi = nut_angles.dpsi; // [rad]
    double deps = nut_angles.deps; // [rad]


    double eps_true = eps_mean + deps;

   
    Matrix Rx_eps_mean = R_x(eps_mean);
    
    
    Matrix Rz_neg_dpsi = R_z(-dpsi);
    
    
    Matrix Rx_neg_eps_true = R_x(-eps_true);

    Matrix NutMat = Rx_neg_eps_true * Rz_neg_dpsi * Rx_eps_mean;
    
    return NutMat;
}