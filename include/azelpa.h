#ifndef AZELPA_H
#define AZELPA_H

#include "matrix.h"
#include "SAT_Const.h"  

/**
 * @file azelpa.h
 * @brief Declaración de la función AzElPa y la estructura de datos para sus resultados.

 */

/**
 * @struct AzElPaData
 * @brief Estructura para almacenar los resultados de la función AzElPa.
 */
struct AzElPaData {
    double Az;      ///< Azimut [rad], en el rango [0, 2*pi)
    double El;      ///< Elevación [rad]
    Matrix dAds;    ///< Parciales de Azimut respecto a s (Matriz 1x3)
    Matrix dEds;    ///< Parciales de Elevación respecto a s (Matriz 1x3)

    // Constructor para inicializar las matrices de parciales a un tamaño conocido
    AzElPaData() : dAds(1, 3), dEds(1, 3) {}
};

/**
 * @brief Calcula azimut, elevación y sus derivadas parciales desde coordenadas locales tangenciales.
 *
 * El sistema de coordenadas local tangencial es Este-Norte-Zenit.
 *
 * @param s Un vector (Matriz 3x1 o 1x3) con las coordenadas topocéntricas locales tangenciales [s_E, s_N, s_Z].
 * Se espera que s(1) sea Este, s(2) Norte, s(3) Zenit.
 * @return AzElPaData Una estructura que contiene Azimut, Elevación y sus parciales.
 **/
AzElPaData AzElPa( Matrix& s);

#endif // AZELPA_H