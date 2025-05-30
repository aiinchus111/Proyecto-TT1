#include "../include/EccAnom.h"

#include <stdexcept>
#include <iostream>

using namespace std;
/**
 * @file eccentric_anomaly.cpp
 * @brief Implementación para calcular la anomalía excéntrica para órbitas elípticas.
 */

/**
 * @brief Calcula la anomalía excéntrica para órbitas elípticas.
 *
 * @param M Anomalía media en [rad]
 * @param e Excentricidad de la órbita [0,1]
 * @return Anomalía excéntrica en [rad]
 * @throws std::runtime_error si la iteración no converge.
 */
double EccAnom(double M, double e) {
    int maxit = 15;
    int i = 1;

    // Valor inicial
    M = fmod(M, 2.0 * Const::pi);
    if (M < 0) {
        M+=2.0*Const::pi;
    }

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = Const::pi;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // Iteración
    while (fabs(f) > 1e-12) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i = i + 1;
        if (i == maxit) {
            cout << ("problemas de convergencia en EccAnom\n") << endl;
            exit(EXIT_FAILURE);
        }
    }

    return E;
}