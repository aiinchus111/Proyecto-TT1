#ifndef ANGLESG_H
#define ANGLESG_H

#include "matrix.h" 

#include "elements.h" 
#include "hgibs.h" 
#include "gibbs.h" 
#include <vector>   

/**
 * @file anglesg.h
 * @brief Declaración de la función anglesg para la determinación inicial de órbita.
 * @author M. Mahooti
 *
 */

/**
 * @struct OrbitStateIOD
 * @brief Estructura para almacenar el resultado de la determinación de órbita (posición y velocidad).
 */
struct OrbitStateIOD {
    Matrix r_vec; ///< Vector de posición (3x1) [m] en ICRF/EME2000 en t2.
    Matrix v_vec; ///< Vector de velocidad (3x1) [m/s] en ICRF/EME2000 en t2.

    OrbitStateIOD() : r_vec(3,1), v_vec(3,1) {} // Inicializa matrices
};

/**
 * @brief Resuelve el problema de determinación de órbita usando tres avistamientos ópticos (método de Gauss).
 *
 * @param az1 Azimut en t1 [rad].
 * @param az2 Azimut en t2 [rad].
 * @param az3 Azimut en t3 [rad].
 * @param el1 Elevación en t1 [rad].
 * @param el2 Elevación en t2 [rad].
 * @param el3 Elevación en t3 [rad].
 * @param Mjd1 Fecha Juliana Modificada de t1 (generalmente UTC para IERS).
 * @param Mjd2 Fecha Juliana Modificada de t2.
 * @param Mjd3 Fecha Juliana Modificada de t3.
 * @param Rs1 Vector de posición (3x1) del sitio de observación 1 en ECEF [m] en Mjd1.
 * @param Rs2 Vector de posición (3x1) del sitio de observación 2 en ECEF [m] en Mjd2.
 * @param Rs3 Vector de posición (3x1) del sitio de observación 3 en ECEF [m] en Mjd3.
 * @return OrbitStateIOD Una estructura que contiene el vector de posición r2 y velocidad v2 en t2.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 * 
 */
OrbitStateIOD anglesg ( double az1, double az2, double az3,
                        double el1, double el2, double el3,
                        double Mjd1, double Mjd2, double Mjd3,
                        Matrix& Rs1_ecef,  Matrix& Rs2_ecef, Matrix& Rs3_ecef );



/** * @brief STUB: Encuentra raíces reales de un polinomio. 
 * @details Coeficientes: poly_coeffs[0]*x^N + poly_coeffs[1]*x^(N-1) + ... + poly_coeffs[N].
 * @param poly_coefficients Vector de coeficientes del polinomio.
 * @return std::vector<double> Vector con las raíces reales encontradas.
 */
std::vector<double> find_real_roots(const std::vector<double>& poly_coefficients);

#endif // ANGLESG_H