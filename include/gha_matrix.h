#ifndef GHA_MATRIX_H
#define GHA_MATRIX_H

#include "matrix.h" 
#include "gast.h"   
#include "r_z.h"    


/**
 * @file gha_matrix.h
 * @brief Declaración de la función para calcular la matriz del Ángulo Horario de Greenwich.
 */

/**
 * @brief Calcula la matriz de transformación desde el ecuador y equinoccio verdaderos
 * al sistema del ecuador terrestre y meridiano de Greenwich.
 *
 * Esta matriz representa la rotación de la Tierra alrededor de su eje polar
 * por un ángulo igual al Tiempo Sidéreo Aparente de Greenwich (GAST).
 *
 * @param Mjd_UT1 Fecha Juliana Modificada en Tiempo Universal UT1.
 * @return Matrix La matriz del Ángulo Horario de Greenwich (3x3).
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif // GHA_MATRIX_H