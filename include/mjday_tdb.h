#ifndef mjday_tdb_H
#define mjday_tdb_H

#include "SAT_Const.h" 

/**
 * @file mjday_tdb.h
 * @brief Utilidades para cálculos de efemérides, como conversiones de tiempo.
 */

/**
 * @brief Convierte Fecha Juliana Modificada en TT a Fecha Juliana Modificada en TDB.
 *
 * Utiliza una serie de aproximación (ej. de Vallado) para la diferencia TDB-TT.
 *
 * @param Mjd_TT Fecha Juliana Modificada en Tiempo Terrestre.
 * @return double Fecha Juliana Modificada en Tiempo Dinámico Baricéntrico.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 * @note Referencia: Vallado D. A; Fundamentals of Astrodynamics and Applications; 3rd ed. (2007).
 */
double Mjday_TDB(double Mjd_TT);


#endif // mjday_tdb_H