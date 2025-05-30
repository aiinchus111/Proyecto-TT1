#ifndef JPL_EPH_DE430_H
#define JPL_EPH_DE430_H

#include "matrix.h" 
#include "Cheb3D.h" 

/**
 * @file jpl_eph_de430.h
 * @brief Calcula posiciones planetarias usando efemérides JPL DE430.

 */

/**
 * @struct PlanetPositions
 * @brief Estructura para almacenar los vectores de posición de los cuerpos celestes.
 * Todos los vectores son 3x1 [x,y,z] en metros, referidos al ICRF.
 * Las posiciones planetarias (Mercurio a Plutón) y Sol son geocéntricas.
 * La posición de la Tierra y la Luna son respecto al SSB (Barycentro del Sistema Solar).
 */
struct PlanetPositions {
    Matrix r_Mercury; Matrix r_Venus; Matrix r_Earth; Matrix r_Mars;
    Matrix r_Jupiter; Matrix r_Saturn; Matrix r_Uranus; Matrix r_Neptune;
    Matrix r_Pluto;   Matrix r_Moon;   Matrix r_Sun;

    // Constructor para inicializar matrices a 3x1
    PlanetPositions() : 
        r_Mercury(3,1), r_Venus(3,1), r_Earth(3,1), r_Mars(3,1),
        r_Jupiter(3,1), r_Saturn(3,1), r_Uranus(3,1), r_Neptune(3,1),
        r_Pluto(3,1),   r_Moon(3,1),   r_Sun(3,1) {}
};

/**
 * @brief Calcula las posiciones ecuatoriales del Sol, la Luna y los nueve planetas principales
 * utilizando las Efemérides JPL DE430.
 *
 * @param Mjd_TDB Fecha Juliana Modificada de TDB (Tiempo Dinámico Baricéntrico).
 * @return PlanetPositions Una estructura con los vectores de posición [m] en el ICRF.
 *
 * @note El tiempo de luz ya está considerado en las efemérides JPL.
 * @note Esta función depende de una variable global `JPL_PC_DATA` (Matrix)
 * que debe ser pre-cargada con los coeficientes de las efemérides.
 * @note También depende de la función 'Cheb3D(...)' para la interpolación de Chebyshev.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2018/01/11   M. Mahooti
 */
PlanetPositions JPL_Eph_DE430(double Mjd_TDB);


#endif // JPL_EPH_DE430_H
