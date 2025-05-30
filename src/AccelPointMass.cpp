#include "../include/AccelPointMass.h"
#include <cmath>   // Para pow() y sqrt()
#include <iostream> // Para la salida de errores

/**
 * @file AccelPointMass.cpp
 * @brief Implementación de la función AccelPointMass.
 */


/**
 * @brief Calcula la aceleración perturbacional sobre un satélite debida a una masa puntual.
 *
 * @param r El vector de posición del satélite (Matriz de 1x3).
 * @param s El vector de posición de la masa puntual (Matriz de 1x3).
 * @param GM El coeficiente gravitacional de la masa puntual.
 * @return Una Matriz que representa el vector de aceleración (Matriz de 1x3).
 */
Matrix& AccelPointMass( Matrix r, Matrix s, double GM) {
    
    if (r.n_column==1) {
        r=r.transpose();
    }
    if (s.n_column==1) {
        s=s.transpose();
    }
    // Verificar las dimensiones de los vectores de entrada
    if (r.n_row != 1 || r.n_column != 3 || s.n_row != 1 || s.n_column != 3) {
        cout << "AccelPointMass: Error - Los vectores de entrada deben ser de 1x3." << endl;
        exit(EXIT_FAILURE);
    }

    // Calcular el vector de posición relativa
    Matrix d = r - s;
    double norm_d = d.norm();
    double norm_s = s.norm();

    // Calcular los términos de la aceleración
    Matrix term1 = d / pow(norm_d, 3);
    Matrix term2 = s / pow(norm_s, 3);

    // Calcular la aceleración
    return (term1 + term2) * (-GM);
}