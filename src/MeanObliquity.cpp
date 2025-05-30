#include "../include/MeanObliquity.h"
#include "../include/SAT_Const.h" 

/**
 * @file mean_obliquity.cpp
 * @brief Implementación de la función para calcular la oblicuidad media de la eclíptica.
 */

double MeanObliquity(double Mjd_TT) {
    // T son siglos Julianos desde J2000.0
    double T = (Mjd_TT - Const::MJD_J2000) / 36525.0;

    double term_T  = 46.8150;
    double term_T2 = 0.00059;
    double term_T3 = 0.001813;

    double epsilon_deg = (84381.448 - (term_T * T + term_T2 * T * T - term_T3 * T * T * T)) / 3600.0;

    double MOblq_rad = epsilon_deg * Const::Rad;

    return MOblq_rad;
}