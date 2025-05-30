#ifndef NUT_MATRIX_H
#define NUT_MATRIX_H

#include "matrix.h"          
#include "meanObliquity.h"   
#include "nutangles.h"        
#include "r_x.h"            
#include "r_z.h"              

/**
 * @file nut_matrix.h
 * @brief Declaración de la función para calcular la matriz de nutación.

 */

/**
 * @brief Calcula la matriz de transformación desde el ecuador y equinoccio medios
 * a los verdaderos (de la fecha).
 *
 * La matriz de nutación se construye como N = Rx(-eps_true) * Rz(-dpsi) * Rx(eps_mean),
 * donde eps_true = eps_mean + deps.
 *
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre).
 * @return Matrix La matriz de nutación 3x3.
 *
 * @note La "ecuación de los equinoccios" dpsi*cos(eps) es la ascensión recta
 * del equinoccio medio referida al ecuador y equinoccio verdaderos y es igual
 * a la diferencia entre el tiempo sidéreo aparente y el medio. 
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix NutMatrix(double Mjd_TT);

#endif // NUT_MATRIX_H