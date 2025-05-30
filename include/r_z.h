#ifndef R_Z_H
#define R_Z_H


#include "matrix.h"

/**
 * @file r_z.h
 * @brief Declaración de la función para generar una matriz de rotación alrededor del eje Z.
 */

/**
 * @brief Genera una matriz de rotación 3x3 alrededor del eje Z.
 *
 * La matriz de rotación se define como:
 * |  cos(a)   sin(a)   0 |
 * | -sin(a)   cos(a)   0 |
 * |    0        0      1 |
 * donde 'a' es el ángulo de rotación.
 *
 * @param angle El ángulo de rotación en radianes.
 * @return Matrix Una matriz 3x3 que representa la rotación.
 */
Matrix R_z(double angle);

#endif // R_Z_H