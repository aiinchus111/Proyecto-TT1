#ifndef cheb3D_H_
#define cheb3D_H_

#include "matrix.h"
#include <stdexcept>

/**
 * @file cheb3D.h
 * @brief Archivo de cabecera para la aproximación de Chebyshev de vectores 3D.
 */

/**
 * @brief Calcula la aproximación de Chebyshev de un vector 3D en un instante de tiempo dado.
 *
 * @param t El instante de tiempo en el que evaluar la aproximación.
 * @param N El número de coeficientes de Chebyshev.
 * @param Ta El inicio del intervalo de tiempo.
 * @param Tb El final del intervalo de tiempo.
 * @param Cx Los coeficientes de Chebyshev para la coordenada x (Matriz 1xN).
 * @param Cy Los coeficientes de Chebyshev para la coordenada y (Matriz 1xN).
 * @param Cz Los coeficientes de Chebyshev para la coordenada z (Matriz 1xN).
 * @return Una Matriz que representa la aproximación del vector 3D en el instante t (Matriz 1x3).
 * @throws std::runtime_error si t está fuera del intervalo [Ta, Tb].
 */
Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz);

#endif // cheb3D_H_