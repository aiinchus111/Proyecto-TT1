#ifndef GEODETIC_H
#define GEODETIC_H

#include "matrix.h" 
#include "SAT_Const.h"  

/**
 * @file geodetic.h
 * @brief Declaración de la función para convertir coordenadas cartesianas a geodésicas.
 * @author M. Mahooti 
 */

/**
 * @struct GeodeticCoords
 * @brief Estructura para almacenar coordenadas geodésicas.
 */
struct GeodeticCoords {
    double lon_rad; ///< Longitud geodésica Este [rad]
    double lat_rad; ///< Latitud geodésica [rad]
    double h_m;     ///< Altitud sobre el elipsoide [m]
};

/**
 * @brief Convierte un vector de posición cartesiano (ECEF) a coordenadas geodésicas.
 *
 * Utiliza un método iterativo para calcular la latitud y la altitud.
 *
 * @param r_ecef Vector de posición cartesiano 3x1 (X, Y, Z) en metros, en el sistema ECEF.
 * @return GeodeticCoords Una estructura que contiene la longitud, latitud y altitud geodésicas.
 * @throws std::runtime_error si el vector de posición r_ecef es el vector cero (norma cero).
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
GeodeticCoords Geodetic(Matrix& r_ecef);

#endif // GEODETIC_H