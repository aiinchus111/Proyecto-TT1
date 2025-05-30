#ifndef GMST_H
#define GMST_H

#include "SAT_Const.h" 
#include "frac.h"  
/**
 * @file gmst.h
 * @brief Declaración de la función para calcular el Tiempo Sidéreo Medio de Greenwich (GMST).
 */

/**
 * @brief Calcula el Tiempo Sidéreo Medio de Greenwich (GMST).
 *
 * Utiliza la fórmula IAU 1982/1994.
 *
 * @param Mjd_UT1 Fecha Juliana Modificada en Tiempo Universal UT1.
 * @return double GMST en radianes, en el rango [0, 2*pi).
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
double gmst(double Mjd_UT1);

#endif // GMST_H