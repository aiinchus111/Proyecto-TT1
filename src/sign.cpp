#include "../include/sign.h"
#include <cmath> 

/**
 * @file sign_util.cpp
 * @brief Implementación de la función sign_.
 */

double sign_(double a, double b) {

    return std::copysign(std::abs(a), b);

}