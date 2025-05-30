#ifndef HGIBBS_H
#define HGIBBS_H

#include "matrix.h" 
#include "angl.h"   
#include "SAT_Const.h"  
#include "gibbs.h" 
#include <string>
#include <vector>

/**
 * @file hgibbs.h
 * @brief Declaración de la función hgibbs para determinación de órbita.
 * @author M. Mahooti (Concepto original MATLAB)
 
 */



/**
 * @brief Implementa el método de Herrick-Gibbs para la determinación de órbitas.
 *
 * Este método determina el vector de velocidad en el punto medio (r2)
 * dados tres vectores de posición (r1, r2, r3) y sus tiempos de observación.
 * Es particularmente útil cuando los ángulos entre los vectores de posición son pequeños.
 *
 * @param r1_in Vector de posición (3x1) #1 en un sistema inercial [m].
 * @param r2_in Vector de posición (3x1) #2 en un sistema inercial [m].
 * @param r3_in Vector de posición (3x1) #3 en un sistema inercial [m].
 * @param Mjd1 Fecha Juliana Modificada de la observación 1.
 * @param Mjd2 Fecha Juliana Modificada de la observación 2.
 * @param Mjd3 Fecha Juliana Modificada de la observación 3.
 * @return GibbsResult Una estructura que contiene v2, theta, theta1, copa y un mensaje de error.
 *
 */
GibbsResult hgibbs ( Matrix& r1_in,  Matrix& r2_in,  Matrix& r3_in,
                    double Mjd1, double Mjd2, double Mjd3);



#endif // HGIBBS_IOD_H