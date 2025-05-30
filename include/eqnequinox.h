#ifndef EQN_EQUINOX_H
#define EQN_EQUINOX_H

#include "nutAngles.h"        
#include "meanObliquity.h"  

/**
 * @file eqn_equinox.h
 * @brief Declaración de la función para calcular la ecuación de los equinoccios.

 */

/**
 * @brief Calcula la ecuación de los equinoccios.
 *
 * La ecuación de los equinoccios (dpsi * cos(eps_mean)) es la ascensión recta
 * del equinoccio medio referida al ecuador y equinoccio verdaderos, y es igual
 * a la diferencia entre el tiempo sidéreo aparente y el medio.
 *
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre).
 * @return double La ecuación de los equinoccios [rad].
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
double EqnEquinox(double Mjd_TT);

#endif // EQN_EQUINOX_H