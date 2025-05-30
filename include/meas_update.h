#ifndef MEAS_UPDATE_H
#define MEAS_UPDATE_H

#include "matrix.h"

/**
 * @file meas_update.h
 * @brief Declaración de la función para la actualización de la medida del filtro de Kalman.
 */

/**
 * @struct MeasUpdateResult
 * @brief Estructura para almacenar los resultados de la actualización de la medida.
 */
struct MeasUpdateResult {
    Matrix K; 
    Matrix x; 
    Matrix P; 

    
    MeasUpdateResult() {} 
    MeasUpdateResult(int n_dim, int m_dim) : K(n_dim, m_dim), x(n_dim, 1), P(n_dim, n_dim) {}
};

/**
 * @brief Realiza la actualización de la medida (paso de corrección) del filtro de Kalman.
 *
 * Calcula la ganancia de Kalman, y actualiza el vector de estado y su matriz de covarianza.
 *
 * @param x_current Vector de estado actual (predicho), x_k^- (n x 1).
 * @param z_measurement Vector de medida actual, z_k (m x 1).
 * @param g_predicted_measurement Medida predicha basada en x_current, h(x_k^-) (m x 1).
 * @param s_std_devs Vector (m x 1 o 1 x m) de desviaciones estándar del ruido de medida.
 * Se usa para construir la matriz diagonal de covarianza del ruido de medida R.
 * @param G_jacobian Matriz de sensibilidad de medida (Jacobiano H_k), (m x n).
 * @param P_current Matriz de covarianza del estado actual (predicha), P_k^- (n x n).
 * @param state_dimension_n Dimensión del vector de estado 'n'.
 * @return MeasUpdateResult Una estructura que contiene K, x actualizado, y P actualizado.
 *
 * @note Concepto original de MATLAB.
 */
MeasUpdateResult MeasUpdate(Matrix& x_current, 
                            Matrix& z_measurement, 
                            Matrix& g_predicted_measurement, 
                            Matrix& s_std_devs, 
                            Matrix& G_jacobian, 
                            Matrix& P_current, 
                            int state_dimension_n);

#endif // MEAS_UPDATE_H