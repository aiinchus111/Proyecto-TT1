#ifndef ACCEL_HARMONIC_H
#define ACCEL_HARMONIC_H

#include "matrix.h"      
#include "legendre.h" 
#include "global.h"

/**
 * @file accel_harmonic.h
 * @brief Declaración de la función para calcular la aceleración debida al campo armónico.

 */

/**
 * @brief Calcula la aceleración debida al campo gravitatorio armónico del cuerpo central.
 *
 * @param r_inertial Vector de posición del satélite (3x1) en el sistema inercial [m].
 * @param E Matriz de transformación (3x3) del sistema inercial al sistema fijo al cuerpo.
 * @param n_max Grado máximo de los armónicos esféricos a considerar.
 * @param m_max Orden máximo de los armónicos esféricos (m_max <= n_max).
 * @param Cnm Matriz (n_max+1 x n_max+1 o similar) con los coeficientes C_nm.
 * Se accede como Cnm(n+1, m+1) para C_{n,m}.
 * @param Snm Matriz (n_max+1 x n_max+1 o similar) con los coeficientes S_nm.
 * Se accede como Snm(n+1, m+1) para S_{n,m}.
 * @return Matrix Vector de aceleración (3x1) en el sistema inercial [m/s^2].
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 * @note Utiliza valores específicos para r_ref (6378.1363e3 m) y gm (398600.4415e9 m^3/s^2)
 * tal como en el script original.
 */
Matrix AccelHarmonic( Matrix& r_inertial,  Matrix& E,
                     int n_max, int m_max);

#endif // ACCEL_HARMONIC_H