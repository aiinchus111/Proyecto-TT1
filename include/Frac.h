#ifndef FRAC_H
#define FRAC_H

#include <cmath> 

/**
 * @file Frac.h
 * @brief Declaración de la función para calcular la parte fraccionaria de un número.
 *
 */

/**
 * @brief Calcula la parte fraccionaria de un número de punto flotante.
 *
 * La función toma un número `x` y devuelve su parte fraccionaria, calculada
 * como `x - std::floor(x)`. El resultado siempre estará en el rango `[0, 1)`.
 * Por ejemplo, Frac(5.75) es 0.75, y Frac(-3.25) es 0.75.
 *
 * @param x El número de punto flotante de entrada.
 * @return double La parte fraccionaria de x. Este valor está siempre en el
 * intervalo semiabierto [0.0, 1.0).
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
double Frac(double x);

#endif // FRAC_H