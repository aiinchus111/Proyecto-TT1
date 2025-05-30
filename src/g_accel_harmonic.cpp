#include "../include/g_accel_harmonic.h"


/**
 * @file g_accel_harmonic.cpp
 * @brief Implementación de la función G_AccelHarmonic.
 */

Matrix G_AccelHarmonic( Matrix& r,  Matrix& U,
                       int n_max, int m_max) {

    const double d_increment = 1.0; // Incremento de posición [m] para diferencias finitas

    Matrix G(3, 3);       // Matriz gradiente, inicializada a ceros si el constructor lo hace
    Matrix dr_perturb(3, 1); // Vector de perturbación

    // Calcular el gradiente por diferencias finitas centrales, columna por columna
    for (int i = 1; i <= 3; ++i) { // Itera para dr_x, dr_y, dr_z (columnas de G)
        
        // Establecer la perturbación dr en la componente i-ésima
        
        // Si Matrix no inicializa a cero por defecto:
        dr_perturb(1,1) = 0.0; dr_perturb(2,1) = 0.0; dr_perturb(3,1) = 0.0;
        dr_perturb(i, 1) = d_increment;

        // Calcular posiciones perturbadas
       
        Matrix r_plus_half_dr  = r + (dr_perturb * (0.5));
        Matrix r_minus_half_dr = r - (dr_perturb * (0.5));

        // Calcular aceleraciones en las posiciones perturbadas
        // AccelHarmonic necesita Cnm y Snm
        Matrix a_plus  = AccelHarmonic(r_plus_half_dr,  U, n_max, m_max);
        Matrix a_minus = AccelHarmonic(r_minus_half_dr, U, n_max, m_max);

        // Diferencia de aceleraciones
        Matrix da_diff = a_plus - a_minus;

        // Derivada con respecto al eje i-ésimo (columna i de G)
        // G(:,i) = da / d_increment
 
        Matrix G_column = da_diff * (1.0 / d_increment);
        
  

        for (int row = 1; row <= 3; ++row) {
            G(row, i) = G_column(row, 1);
        }
    }

    return G;
}