#ifndef DEINTEG_H
#define DEINTEG_H

#include "matrix.h" 
#include <vector>
#include <string>
#include <functional>

/**
 * @file deinteg.h
 * @brief Integrador de Ecuaciones Diferenciales Ordinarias (método de Shampine & Gordon).
 */


enum class DE_Status {
    DE_INIT = 1,      // Reiniciar integración / Primera llamada
    DE_DONE = 2,      // Paso/Integración completada exitosamente
    DE_BADACC = 3,    // No se pudo alcanzar la precisión requerida
    DE_NUMSTEPS = 4,  // Número máximo de pasos excedido
    DE_STIFF = 5,     // Se sospecha problema rígido (stiff)
    DE_INVPARAM = 6   // Parámetros de entrada inválidos
};

// Estructura para devolver el resultado 
struct DEIntegResult {
    Matrix y_out;        // Vector de estado en tout
    double t_out;        // Tiempo final alcanzado (podría ser < tout si hay error)
    DE_Status status;    // Estado final de la integración
    

    DEIntegResult(int n_eqn) : y_out(n_eqn, 1), t_out(0.0), status(DE_Status::DE_INIT) {}
};


/**
 * @brief Integra un sistema de ecuaciones diferenciales ordinarias dy/dt = func(t,y).
 * Implementa el método multipaso de orden variable y tamaño de paso variable de Shampine & Gordon.
 *
 * @param func Función que calcula las derivadas dy/dt.
 * Debe tener la firma: Matrix func(double t, const Matrix& y_current)
 * @param t Tiempo inicial de la integración.
 * @param tout Tiempo final deseado para la integración.
 * @param relerr Tolerancia de error relativo.
 * @param abserr Tolerancia de error absoluto.
 * @param n_eqn Número de ecuaciones en el sistema (dimensión de y).
 * @param y_in Vector de estado inicial (Matrix n_eqn x 1).
 * @return Matrix Vector de estado y en el tiempo tout (o en el último tiempo alcanzado si hay error).
 
 *

 * @note Basado en Shampine, Gordon: "Computer solution of Ordinary Differential Equations".
 * @note Last modified (MATLAB): 2015/08/25 M. Mahooti
 */
Matrix DEInteg(std::function<Matrix(double, const Matrix&)> func,
               double t, double tout,
               double relerr, double abserr,
               int n_eqn, Matrix y_in); 

#endif // DEINTEG_H