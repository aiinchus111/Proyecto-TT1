#ifndef IERS_H
#define IERS_H

#include "matrix.h" 
#include "SAT_Const.h"  

/**
 * @file iers.h
 * @brief Definición de la estructura EopResults y declaración de la función IERS.
 */

/**
 * @struct EopResults
 * @brief Estructura para almacenar los Parámetros de Orientación Terrestre (EOP).
 */
struct EopResults {
    double x_pole;    
    double y_pole;   
    double UT1_UTC;  
    double LOD;       
    double dpsi;     
    double deps;      
    double dx_pole;   
    double dy_pole;   
    double TAI_UTC;   
};

/**
 * @brief Gestiona y recupera los datos de tiempo y movimiento polar del IERS.
 *
 * La matriz `eop_data` debe contener los datos EOP con la siguiente estructura de filas (1-indexada):
 * - Fila 4: MJD (entero o con fracción de día)
 * - Fila 5: x_pole (en arcosegundos)
 * - ...
 * - Fila 13: TAI-UTC (en segundos)
 *
 * @param eop_matrix Matriz  que contiene los datos EOP.
 * @param Mjd_UTC Fecha Juliana Modificada en UTC para la cual se requieren los parámetros.
 * @param interp_option Carácter que indica el tipo de interpolación:
 * 'l' para interpolación lineal.
 * 'n' para usar los valores del día más cercano.
 * Valor por defecto es 'n'.
 * @return EopResults Una estructura con los parámetros EOP calculados/interpolados.
 * @throws std::runtime_error Si el MJD no se encuentra o datos insuficientes para interpolar.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2018/02/01   M. Mahooti
 */
EopResults IERS( Matrix& eop_matrix, double Mjd_UTC, char interp_option = 'n');

#endif // IERS_H