#ifndef VAR_EQN_H
#define VAR_EQN_H

#include "matrix.h"      
#include "iers.h"
#include "timediff.h"
#include "precc_matrix.h"
#include "nut_matrix.h"
#include "polematrix.h"
#include "gha_matrix.h"
#include "accelharmonic.h"
#include "g_accel_harmonic.h" // Para el gradiente
#include "SAT_Const.h"            // Para constantes globales Const::
#include "global.h"
/**
 * @file var_eqn.h
 * @brief Declaración de la función para calcular las ecuaciones variacionales.
 * @author M. Mahooti (Concepto original MATLAB)
 * @author Tu Nombre/Alias (Adaptación C++)
 * @version 1.0
 * @date 2025-05-22
 */

/**
 * @brief Calcula las ecuaciones variacionales: la derivada del vector de estado
 * y de la matriz de transición de estados (STM).
 *
 * @param time_from_epoch_sec Tiempo en segundos desde la época de referencia definida en aux_params.
 * @param yPhi_augmented_state Vector columna (42x1) que contiene el vector de estado
 * (6 elementos: rx,ry,rz,vx,vy,vz) seguido de la matriz de transición de estados Phi (36 elementos,
 * almacenada por columnas).
 * @param aux_params Estructura con parámetros auxiliares de la simulación, incluyendo
 * épocas de referencia MJD_UTC_epoch, Mjd_TT_epoch, parámetros del modelo de gravedad (n, m, Cnm, Snm).
 * @return Matrix Derivada del vector yPhi_augmented_state (42x1).
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 * @note Esta función depende de variables globales pre-cargadas para EOP (`eopdata`).
 * @note Los coeficientes Cnm y Snm se esperan dentro de aux_params.
 */
Matrix VarEqn(double time_from_epoch_sec, 
              Matrix& yPhi_augmented_state);

#endif // VAR_EQN_H