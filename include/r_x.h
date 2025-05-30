#ifndef R_X_H
#define R_X_H

#include "matrix.h"

/**
 * @file r_x.h
 * @brief Declaración de la función para generar una matriz de rotación alrededor del eje X.

 */

/**
 * @brief Genera una matriz de rotación 3x3 alrededor del eje X.
 *
 * La matriz de rotación se define como:
 * | 1    0        0      |
 * | 0  cos(a)   sin(a) |
 * | 0 -sin(a)   cos(a) |
 * donde 'a' es el ángulo de rotación.
 *
 * @param angle El ángulo de rotación en radianes.
 * @return Matrix Una matriz 3x3 que representa la rotación.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix R_x(double angle);

#endif 