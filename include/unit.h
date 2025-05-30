#ifndef UNIT_H
#define UNIT_H

#include "matrix.h" 
#include <cmath>   

/**
 * @file unit.h
 * @brief Declaración de la función para calcular un vector unitario.
 * @author M. Mahooti (Concepto original MATLAB)
 
 */

/**
 * @brief Calcula un vector unitario a partir del vector original.
 *
 * Si la magnitud del vector de entrada es menor que un umbral pequeño,
 * se devuelve un vector cero. De lo contrario, se devuelve el vector
 * normalizado (unitario).
 *
 * @param vec El vector de entrada (se asume un objeto Matrix, típicamente 3x1).
 * @return Matrix El vector unitario (3x1) o un vector cero (3x1).
 *

 */
Matrix unit( Matrix& vec);

#endif // UNIT_VECTOR_H