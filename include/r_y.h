#ifndef R_Y_H
#define R_Y_H

#include "matrix.h" 

/**
 * @file r_y.h
 * @brief Declaración de la función para generar una matriz de rotación alrededor del eje Y.
 */

/**
 * @brief Genera una matriz de rotación 3x3 alrededor del eje Y.
 *
 * La matriz de rotación se define como:
 * |  cos(a)   0   -sin(a) |
 * |    0      1      0    |
 * |  sin(a)   0    cos(a) |
 * donde 'a' es el ángulo de rotación.
 *
 * @param angle El ángulo de rotación en radianes.
 * @return Matrix Una matriz 3x3 que representa la rotación.
 *
 */
Matrix R_y(double angle);

#endif // R_Y_H