#ifndef G_ACCEL_HARMONIC_H
#define G_ACCEL_HARMONIC_H

#include "matrix.h"          
#include "accelharmonic.h"  

/**
 * @file g_accel_harmonic.h
 * @brief Declaración de la función para calcular el gradiente del campo gravitatorio armónico.
 */

/**
 * @brief Calcula el gradiente del campo gravitatorio armónico de la Tierra (o cuerpo central).
 *
 * El gradiente (G = da/dr) se calcula en el sistema verdadero de la fecha (true-of-date)
 * utilizando un método de diferencias finitas centrales.
 *
 * @param r Vector de posición del satélite (3x1) en el sistema verdadero de la fecha [m].
 * @param U Matriz de transformación (3x3) desde el sistema verdadero de la fecha 
 * al sistema fijo al cuerpo. (Equivalente a 'E' en AccelHarmonic).
 * @param n_max Grado máximo del modelo gravitatorio.
 * @param m_max Orden máximo del modelo gravitatorio (m_max <= n_max).
 * @param Cnm Matriz que contiene los coeficientes C_nm del modelo gravitatorio.
 * @param Snm Matriz que contiene los coeficientes S_nm del modelo gravitatorio.
 * @return Matrix La matriz gradiente G (3x3) [1/s^2].
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix G_AccelHarmonic( Matrix& r,  Matrix& U,
                       int n_max, int m_max);

#endif // G_ACCEL_HARMONIC_H