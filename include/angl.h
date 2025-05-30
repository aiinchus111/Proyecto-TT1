#ifndef ANGL_H
#define ANGL_H

#include "matrix.h" 
#include <cmath>   

/**
 * @file angl.h
 * @brief Declaración de la función para calcular el ángulo entre dos vectores.
 * @author M. Mahooti 
 */

/**
 * @brief Calcula el ángulo (en radianes) entre dos vectores.
 *
 * El ángulo devuelto está en el rango [0, pi].
 *
 * @param vec1 Primer vector (objeto Matrix, se asume 3x1 o similar).
 * @param vec2 Segundo vector (objeto Matrix, se asume 3x1 o similar y de las mismas dimensiones que vec1).
 * @return double El ángulo theta entre los dos vectores en radianes.
 * Devuelve un valor numérico grande (999999.1) si el ángulo es indefinido
 */
double angl( Matrix& vec1,  Matrix& vec2);



#endif // ANGL_H