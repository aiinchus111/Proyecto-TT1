#ifndef GAST_H
#define GAST_H

#include "gmst.h"          
#include "eqnequinox.h"   
#include "SAT_Const.h"        

/**
 * @file gast.h
 * @brief Declaración de la función para calcular el Tiempo Sidéreo Aparente de Greenwich (GAST).

 */

/**
 * @brief Calcula el Tiempo Sidéreo Aparente de Greenwich (GAST).
 *
 * GAST = GMST + Ecuación de los Equinoccios, normalizado a [0, 2*pi).
 *
 * @param Mjd_UT1 Fecha Juliana Modificada en Tiempo Universal UT1.
 * @return double GAST en radianes, en el rango [0, 2*pi).
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
double gast(double Mjd_UT1);

#endif // GAST_H