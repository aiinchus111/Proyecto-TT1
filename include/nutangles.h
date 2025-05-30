#ifndef NUTANGLES_H
#define NUTANGLES_H

#include "SAT_Const.h"

/**
 * @file nutation.h
 * @brief Declaración de la función NutAngles y la estructura para sus resultados.
 * @author M. Mahooti 
 */

/**
 * @struct NutationAnglesResult
 * @brief Estructura para almacenar los ángulos de nutación calculados.
 */
struct NutationAnglesResult {
    double dpsi; ///< Nutación en longitud [rad]
    double deps; ///< Nutación en oblicuidad [rad]
};

/**
 * @brief Calcula los ángulos de nutación en longitud (dpsi) y oblicuidad (deps).
 *
 * Utiliza la serie de nutación IAU 1980.
 *
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre).
 * @return NutationAnglesResult Una estructura que contiene dpsi y deps en radianes.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
NutationAnglesResult NutAngles(double Mjd_TT);

#endif // NUTANGLES_H