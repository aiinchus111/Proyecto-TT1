#ifndef PRECC_MATRIX_H
#define PRECC_MATRIX_H

#include "matrix.h" 
#include "SAT_Const.h"  
#include "r_y.h"    
#include "r_z.h"    

/**
 * @file precc_matrix.h
 * @brief Declaración de la función para calcular la matriz de precesión IAU 1976.
 */

/**
 * @brief Calcula la matriz de transformación de precesión para coordenadas ecuatoriales.
 *
 * Transforma coordenadas desde la época Mjd_1 a la época Mjd_2 utilizando el
 * modelo de precesión IAU 1976 (Lieske et al., 1977).
 *
 * @param Mjd_1 Época inicial (Fecha Juliana Modificada, Tiempo Terrestre).
 * @param Mjd_2 Época a la que se precesiona (Fecha Juliana Modificada, Tiempo Terrestre).
 * @return Matrix La matriz de precesión 3x3.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif // PRECC_MATRIX_H