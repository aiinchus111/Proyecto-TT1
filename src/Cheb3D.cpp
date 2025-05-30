#include "../include/Cheb3D.h"
#include <stdexcept> 

/**
 * @file cheb3D.cpp
 * @brief Implementación de la función Cheb3D.
 */

/**
 * @brief Calcula la aproximación de Chebyshev de un vector 3D en un instante de tiempo dado.
 *
 * @param t El instante de tiempo en el que evaluar la aproximación.
 * @param N El número de coeficientes de Chebyshev.
 * @param Ta El inicio del intervalo de tiempo.
 * @param Tb El final del intervalo de tiempo.
 * @param Cx Los coeficientes de Chebyshev para la coordenada x (Matriz 1xN).
 * @param Cy Los coeficientes de Chebyshev para la coordenada y (Matriz 1xN).
 * @param Cz Los coeficientes de Chebyshev para la coordenada z (Matriz 1xN).
 * @return Una Matriz que representa la aproximación del vector 3D en el instante t (Matriz 1x3).
 */
Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz) {
    // Comprobar validez
    if ((t < Ta) || (Tb < t)) {
        cout<<("ERROR: Tiempo fuera de rango en Cheb3D::Value\n");
        exit(EXIT_FAILURE);
    }

    // Algoritmo de Clenshaw
    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1(1, 3);
    Matrix f2(1, 3);
    Matrix aux(3);

    for (int i = N; i >= 2; --i) {
        Matrix old_f1 = f1;
        aux(1)=Cx(i); aux(2)=Cy(i); aux(3)=Cz(i);
        f1 =  f1 *2 * tau - f2 + aux;
        f2 = old_f1;
    }
    aux(1)=Cx(1); aux(2)=Cy(1); aux(3)=Cz(1);
    return   f1 *tau - f2 + aux;
}