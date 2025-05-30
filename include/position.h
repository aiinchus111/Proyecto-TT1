#ifndef POSITION_H
#define POSITION_H

#include "matrix.h"
#include "SAT_Const.h" 
/**
 * @file position.h
 * @brief Declaración de la función para convertir coordenadas geodésicas a cartesianas.
 */

/**
 * @brief Calcula el vector de posición cartesiano (ECEF) a partir de coordenadas geodésicas.
 *
 * Convierte longitud, latitud geodésica y altitud sobre el elipsoide
 * a un vector de posición cartesiano tridimensional [x, y, z].
 *
 * @param lon Longitud geodésica en radianes.
 * @param lat Latitud geodésica en radianes.
 * @param h Altitud sobre el elipsoide de referencia en metros.
 * @return Matrix Un vector columna (3x1) con las coordenadas cartesianas [x, y, z] en metros.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix Position(double lon, double lat, double h);

#endif // POSITION_H