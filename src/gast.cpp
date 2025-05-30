#include "../include/gast.h"
#include <cmath> 

/**
 * @file gast.cpp
 * @brief Implementación de la función para calcular el Tiempo Sidéreo Aparente de Greenwich (GAST).
 */

// Función auxiliar para el módulo, replicando el comportamiento de mod(a,n) de MATLAB para n > 0
// Devuelve un resultado en [0, divisor)
// Esta función podría estar en un archivo de utilidades comunes si se usa en múltiples sitios.
static double custom_mod_positive(double value, double divisor) {
    if (divisor == 0.0) {
    
        return value; 
    }
    double result = std::fmod(value, divisor);
    if (result < 0.0) {
        result += divisor;
    }
    return result;
}

double gast(double Mjd_UT1) {
    // Calcular GMST
    double gmst_val = gmst(Mjd_UT1); // [rad]

    // Calcular la Ecuación de los Equinoccios
    double eqn_eq_val = EqnEquinox(Mjd_UT1); // [rad]

    // Calcular GAST
    double gast_val = gmst_val + eqn_eq_val;

    // Normalizar a [0, 2*pi)
    // Const::pi2 debe estar definido en const.h como 2.0 * Const::pi
    double gast_normalized = custom_mod_positive(gast_val, Const::pi2);
    
    return gast_normalized;
}