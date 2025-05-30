#ifndef GIBBS_IOD_H
#define GIBBS_IOD_H

#include "matrix.h" 
#include "angl.h"   
#include "SAT_Const.h"  
#include <string>
#include <vector>  

/**
 * @file gibbs.h
 * @brief Declaración de la función gibbs para determinación de órbita.
 * @author M. Mahooti (Concepto original MATLAB)
 */

/**
 * @struct GibbsResult
 * @brief Estructura para almacenar los resultados del método de Gibbs.
 */
struct GibbsResult {
    Matrix v2;          ///< Vector de velocidad (3x1) en r2 [m/s].
    double theta;       ///< Ángulo entre r1 y r2 [rad].
    double theta1;      ///< Ángulo entre r2 y r3 [rad].
    double copa;        ///< Ángulo de coplanaridad (entre r1 y el plano r2-r3) [rad].
    std::string error_msg; ///< Mensaje de estado: "ok", "not coplanar", "impossible".

    /** @brief Constructor por defecto, inicializa v2 a 3x1 ceros y error a "ok". */
    GibbsResult() : v2(3,1), theta(0.0), theta1(0.0), copa(0.0), error_msg("           ok") {}
};

/**
 * @brief Implementa el método de Gibbs para la determinación de órbitas.
 * Determina el vector de velocidad en el punto medio (r2) dados tres
 * vectores de posición coplanares r1, r2, r3.
 *
 * @param r1 Vector de posición (3x1) #1 en un sistema inercial [m].
 * @param r2 Vector de posición (3x1) #2 en un sistema inercial [m].
 * @param r3 Vector de posición (3x1) #3 en un sistema inercial [m].
 * @return GibbsResult Una estructura que contiene v2, theta, theta1, copa y un mensaje de error.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti

 */
GibbsResult gibbs( Matrix& r1,  Matrix& r2,  Matrix& r3);


Matrix unit_vector_helper( Matrix& vec);


#endif // GIBBS_IOD_H