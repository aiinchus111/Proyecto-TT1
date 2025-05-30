#include "../include/precc_matrix.h"
#include <cmath> 
/**
 * @file precC_matrix.cpp
 * @brief Implementación de la función para calcular la matriz de precesión IAU 1976.
 */

Matrix PrecMatrix(double Mjd_1, double Mjd_2) {
    // Tiempo T en siglos Julianos desde J2000.0 hasta la época inicial Mjd_1
    double T  = (Mjd_1 - Const::MJD_J2000) / 36525.0;
    // Intervalo de tiempo dT entre Mjd_1 y Mjd_2 en siglos Julianos
    double dT = (Mjd_2 - Mjd_1) / 36525.0;

    // Polinomios para los ángulos de precesión (coeficientes en arcosegundos)
    // Zeta_A
    double P0_T_zeta = 2306.2181 + (1.39656 - 0.000139 * T) * T;
    double P1_T_zeta = 0.30188 - 0.000344 * T;
    double P2_zeta   = 0.017998;
    double zeta_arcsec = (P0_T_zeta + (P1_T_zeta + P2_zeta * dT) * dT) * dT;
    double zeta_rad = zeta_arcsec / Const::Arcs;

    // z_A
    double S0_T_z = 0.79280 + 0.000411 * T;
    double S1_z   = 0.000205;
    double z_arcsec_offset = (S0_T_z + S1_z * dT) * dT * dT;
    double z_rad = zeta_rad + (z_arcsec_offset / Const::Arcs);

    // Theta_A
    double Q0_T_theta = 2004.3109 - (0.85330 + 0.000217 * T) * T;
    double Q1_T_theta = 0.42665 + 0.000217 * T;
    double Q2_theta   = 0.041833;
    double theta_arcsec = (Q0_T_theta - (Q1_T_theta + Q2_theta * dT) * dT) * dT;
    double theta_rad = theta_arcsec / Const::Arcs;

    // Construir las matrices de rotación elemental
    Matrix Rz_neg_z     = R_z(-z_rad);
    Matrix Ry_theta     = R_y(theta_rad);
    Matrix Rz_neg_zeta  = R_z(-zeta_rad);

    // Calcular la matriz de precesión
    // PrecMat = R_z(-z) * R_y(theta) * R_z(-zeta)
    Matrix PrecMat = Rz_neg_z * Ry_theta * Rz_neg_zeta;
    
    return PrecMat;
}