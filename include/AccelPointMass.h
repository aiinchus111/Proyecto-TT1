#ifndef POINT_MASS_ACCELERATION_H_
#define POINT_MASS_ACCELERATION_H_

#include "matrix.h"

/**
 * @file point_mass_acceleration.h
 * @brief Archivo de cabecera para el cálculo de la aceleración perturbacional debida a una masa puntual.
 */

/**
 * @brief Calcula la aceleración perturbacional sobre un satélite debida a una masa puntual.
 *
 * @param r El vector de posición del satélite (Matriz de 1x3).
 * @param s El vector de posición de la masa puntual (Matriz de 1x3).
 * @param GM El coeficiente gravitacional de la masa puntual.
 * @return Una Matriz que representa el vector de aceleración (Matriz de 1x3).
 */
Matrix& AccelPointMass(Matrix r, Matrix s, double GM);

#endif // POINT_MASS_ACCELERATION_H_