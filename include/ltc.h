#ifndef LTC_H
#define LTC_H

#include "matrix.h" 
#include "r_y.h"   
#include "r_z.h"    

/**
 * @file ltc.h
 * @brief Declaración de la función LTC para la transformación a coordenadas locales tangenciales.
 */

/**
 * @brief Calcula la matriz de transformación desde el sistema del meridiano de Greenwich (ECEF)
 * a coordenadas locales tangenciales (ENZ - East, North, Zenith).
 *
 * La transformación se realiza mediante una secuencia de rotaciones y una permutación de filas.
 * La matriz resultante M transforma un vector de ECEF a ENZ: v_ENZ = M * v_ECEF.
 *
 * @param lon Longitud geodésica Este [rad].
 * @param lat Latitud geodésica [rad].
 * @return Matrix Matriz de rotación 3x3.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix LTC(double lon, double lat);

#endif // LTC_H