#ifndef POLEMATRIX_H
#define POLEMATRIX_H

#include "matrix.h" 
#include "r_x.h"    
#include "r_y.h"    

/**
 * @file polematrix.h
 * @brief Declaración de la función para calcular la matriz de movimiento polar.
 */

/**
 * @brief Calcula la matriz de transformación desde coordenadas pseudo-fijas a la Tierra
 * a coordenadas fijas a la Tierra para una fecha dada, debido al movimiento del polo.
 *
 * La matriz de movimiento polar W se calcula como W = Ry(-xp) * Rx(-yp).
 *
 * @param xp Coordenada x del polo (CIP en el sistema terrestre) en radianes.
 * @param yp Coordenada y del polo (CIP en el sistema terrestre) en radianes.
 * @return Matrix La matriz de movimiento polar 3x3.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix PoleMatrix(double xp, double yp);

#endif // POLEMATRIX_H