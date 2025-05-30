#ifndef ECCANOM
#define ECCANOM

#include <cmath>
#include <stdexcept>
#include "SAT_Const.h"

/**
 * @file EccAnom.h
 * @brief Archivo de cabecera para calcular la anomalía excéntrica para órbitas elípticas.
 */

/**
 * @brief Calcula la anomalía excéntrica para órbitas elípticas.
 *
 * @param M Anomalía media en [rad]
 * @param e Excentricidad de la órbita [0,1]
 * @return Anomalía excéntrica en [rad]
 */
double EccAnom(double M, double e);

#endif // ECCANOM