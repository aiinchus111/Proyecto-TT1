#ifndef ORBITAL_ELEMENTS_H
#define ORBITAL_ELEMENTS_H

#include "matrix.h" 
#include "SAT_Const.h"  

/**
 * @file elements.h
 * @brief Declaración de la función para calcular elementos orbitales Keplerianos.
 * @author M. Mahooti 

 */

/**
 * @struct OrbitalElementsResult
 * @brief Estructura para almacenar los elementos orbitales Keplerianos osculadores.
 */
struct OrbitalElements {
    double p_semilatus_rectum; ///< Semilatus rectum [m]
    double a_semimajor_axis;   ///< Semieje mayor [m]
    double e_eccentricity;     ///< Excentricidad (adimensional)
    double i_inclination_rad;  ///< Inclinación [rad], rango [0, pi]
    double Omega_raan_rad;     ///< Longitud del nodo ascendente (RAAN) [rad], rango [0, 2*pi)
    double omega_argp_rad;     ///< Argumento del pericentro [rad], rango [0, 2*pi)
    double M_mean_anomaly_rad; ///< Anomalía media [rad], rango [0, 2*pi)

    // Se podrían añadir otros elementos si se desea, como nu (anomalía verdadera), E (anomalía excéntrica)
    double nu_true_anomaly_rad;///< Anomalía verdadera [rad]
    double E_ecc_anomaly_rad;  ///< Anomalía excéntrica [rad]
    double u_arg_latitude_rad; ///< Argumento de la latitud [rad]


    /** @brief Constructor por defecto, inicializa a cero. */
    OrbitalElements () : p_semilatus_rectum(0), a_semimajor_axis(0), e_eccentricity(0),
                              i_inclination_rad(0), Omega_raan_rad(0), omega_argp_rad(0),
                              M_mean_anomaly_rad(0), nu_true_anomaly_rad(0),
                              E_ecc_anomaly_rad(0), u_arg_latitude_rad(0) {}
};

/**
 * @brief Calcula los elementos orbitales Keplerianos osculadores a partir del vector de estado.
 * Válido para órbitas elípticas.
 *
 * @param y_state Vector de estado 6x1 [rx, ry, rz, vx, vy, vz]^T en un sistema inercial (ej. ICRF/EME2000).
 * Las unidades deben ser metros para posición y m/s para velocidad.
 * @return OrbitalElementsResult Una estructura que contiene los 7 elementos orbitales principales.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 */
OrbitalElements elements( Matrix& y_state);


#endif // ORBITAL_ELEMENTS_H