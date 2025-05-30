#ifndef TIME_UPDATE_H
#define TIME_UPDATE_H

#include "matrix.h" 

/**
 * @file time_update.h
 * @brief Declaración de las funciones para la actualización temporal de la matriz de covarianza.
 */

/**
 * @brief Realiza la actualización temporal de la matriz de covarianza (propagación).
 *
 * Calcula P_k^- = Phi * P_{k-1}^+ * Phi^T + Q_{k-1}.
 *
 * @param P_current La matriz de covarianza actual P_{k-1}^+ (entrada).
 * @param Phi La matriz de transición de estados Phi_{k-1}.
 * @param Qdt La matriz de covarianza del ruido del proceso Q_{k-1}.
 * @return Matrix La matriz de covarianza propagada P_k^-.
 *
 */
Matrix TimeUpdate( Matrix& P_current,  Matrix& Phi,  Matrix& Qdt);

/**
 * @brief Sobrecarga de TimeUpdate sin matriz de covarianza del ruido del proceso Qdt.
 *
 * Calcula P_k^- = Phi * P_{k-1}^+ * Phi^T.
 *
 * @param P_current La matriz de covarianza actual P_{k-1}^+ (entrada).
 * @param Phi La matriz de transición de estados Phi_{k-1}.
 * @return Matrix La matriz de covarianza propagada P_k^-.
 */
Matrix TimeUpdate(Matrix& P_current,  Matrix& Phi);

#endif // TIME_UPDATE_H