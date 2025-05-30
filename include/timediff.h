#ifndef TIMEDIFF_H
#define TIMEDIFF_H

/**
 * @file timediff.h
 * @brief Definición de la estructura para los resultados de timediff y declaración de la función.
 */

/**
 * @struct TimeDifferences
 * @brief Estructura para almacenar las diferencias de tiempo calculadas por la función timediff.
 * Todas las diferencias están en segundos.
 */
struct TimeDifferences {
    double UT1_TAI; ///< Diferencia UT1 - TAI [s]
    double UTC_GPS; ///< Diferencia UTC - GPS Time [s]
    double UT1_GPS; ///< Diferencia UT1 - GPS Time [s]
    double TT_UTC;  ///< Diferencia TT - UTC [s]
    double GPS_UTC; ///< Diferencia GPS Time - UTC [s]
};

/**
 * @brief Calcula varias diferencias entre escalas de tiempo.
 *
 * Las escalas de tiempo involucradas son UT1 (Tiempo Universal), UTC (Tiempo Universal Coordinado),
 * TAI (Tiempo Atómico Internacional), TT (Tiempo Terrestre) y GPS Time.
 *
 * @param UT1_UTC Diferencia UT1 - UTC en segundos.
 * @param TAI_UTC Diferencia TAI - UTC en segundos (segundos intercalares).
 * @return TimeDifferences Una estructura que contiene las diferencias de tiempo calculadas.
 *
 * @note Concepto original de MATLAB:
 * Last modified:   2015/08/12   M. Mahooti
 * @note Constantes internas usadas:
 * TT_TAI  = +32.184 s
 * GPS_TAI = -19.0   s
 */
TimeDifferences timediff(double UT1_UTC, double TAI_UTC);

#endif // TIMEDIFF_H