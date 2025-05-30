#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "matrix.h" 

/**
 * @file legendre.h
 * @brief Declaración de la función para calcular polinomios de Legendre asociados y sus derivadas.

 */

/**
 * @struct LegendreOutput
 * @brief Estructura para almacenar los resultados de la función Legendre.
 * Contiene las matrices de los polinomios de Legendre asociados (pnm)
 * y sus derivadas respecto a fi (dpnm).
 */
struct LegendreOutput {
    Matrix pnm;  
    Matrix dpnm; 

    /**
     * @brief Constructor para LegendreOutput.
     * Inicializa las matrices pnm y dpnm con las dimensiones dadas.
     * @param n_max Grado máximo.
     * @param m_max Orden máximo.
     */
    LegendreOutput(int n_max, int m_max) : pnm(n_max + 1, m_max + 1), dpnm(n_max + 1, m_max + 1) {
    }
};

/**
 * @brief Calcula los polinomios asociados de Legendre normalizados P_nm(sin(fi))
 * y sus primeras derivadas dP_nm(sin(fi))/dfi.
 *
 * Las matrices de salida `pnm` y `dpnm` son de tamaño (n+1)x(m+1).
 * El elemento `pnm(i+1, j+1)` corresponde a P_ij (grado i, orden j).
 * Similarmente para `dpnm`.
 *
 * @param n Grado máximo de los polinomios.
 * @param m Orden máximo de los polinomios (m <= n).
 * @param fi Argumento angular phi en radianes. Los polinomios son funciones de sin(fi).
 * @return LegendreOutput Una estructura que contiene las matrices pnm y dpnm.
 * @throws std::invalid_argument si m > n.
 *
 * @note Concepto original de MATLAB. Las fórmulas corresponden a polinomios
 * asociados de Legendre completamente normalizados, comunes en geofísica.
 */
LegendreOutput Legendre(int n, int m, double fi);

#endif // LEGENDRE_H