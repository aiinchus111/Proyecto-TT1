#ifndef ACCEL_H
#define ACCEL_H

#include "matrix.h"       


/**
 * @file accel.h
 * @brief Declaración de la función principal de cálculo de aceleración.
 
 */

/**
 * @brief Calcula la aceleración total de un satélite en órbita terrestre.
 *
 * Incluye el campo gravitatorio armónico de la Tierra, perturbaciones
 * gravitatorias del Sol y la Luna, y opcionalmente otros planetas.
 * El resultado es la derivada del vector de estado dY = [velocidad; aceleracion].
 *
 * @param time_from_epoch_sec Tiempo en segundos desde la época de referencia definida en aux_params.Mjd_UTC_epoch.
 * @param Y_state Vector de estado del satélite (6x1) [rx,ry,rz,vx,vy,vz]^T en el sistema ICRF/EME2000 [m, m/s].
 * @param aux_params Estructura con parámetros auxiliares de la simulación (época, flags, modelo de gravedad).
 * @return Matrix Derivada del vector de estado dY (6x1) en el sistema ICRF/EME2000 [m/s; m/s^2].
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti (MATLAB original; última modif. de C++ es la fecha actual)
 * @note Esta función depende de variables globales pre-cargadas para EOP (`eopdata`)
 * y para coeficientes de efemérides JPL (`JPL_PC_DATA`).
 * @note Los coeficientes Cnm y Snm para AccelHarmonic se esperan dentro de aux_params.
 */
Matrix Accel(double time_from_epoch_sec, 
             Matrix& Y_state);
             

#endif // ACCEL_H