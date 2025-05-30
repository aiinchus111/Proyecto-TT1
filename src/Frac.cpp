#include "../include/Frac.h" 
#include <cmath>  

/**
 * @file Frac.cpp
 * @brief Implementación de la función para calcular la parte fraccionaria de un número.
 */


double Frac(double x) {
    return x - std::floor(x);
}