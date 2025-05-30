#ifndef MEANOBLIQUITY_H
#define MEANOBLIQUITY_H

/**
 * @file MeanObliquity.h
 * @brief Declaración de la función para calcular la oblicuidad media de la eclíptica.
 */

/**
 * @brief Calcula la oblicuidad media de la eclíptica.
 *
 * Basado en la fórmula de Laskar (1986) tal como se presenta en el
 * "Explanatory Supplement to the Astronomical Almanac" (1992), p. 703, fórmula 11.22.
 *
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre).
 * @return double La oblicuidad media de la eclíptica en radianes.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
double MeanObliquity(double Mjd_TT);

#endif // MEANOBLIQUITY_H