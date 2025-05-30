#include "../include/unit.h"
#include <stdexcept> 

/**
 * @file unit_vector.cpp
 * @brief Implementación de la función para calcular un vector unitario.
 */

Matrix unit( Matrix& vec) {
    
   
    const double small_threshold = 1.0e-6; // Umbral como en el código MATLAB

    Matrix outvec(3, 1); // Inicializar el vector de salida (asumiendo que se inicializa a ceros)

    double magv = vec.norm(); // Asume que vec.norm() existe y calcula la norma L2

    if (magv > small_threshold) {
        
        double inv_magv = 1.0 / magv;
        outvec(1,1) = vec(1,1) * inv_magv;
        outvec(2,1) = vec(2,1) * inv_magv;
        outvec(3,1) = vec(3,1) * inv_magv;
    } else {
      
        outvec(1,1) = 0.0;
        outvec(2,1) = 0.0;
        outvec(3,1) = 0.0;
    }

    return outvec;
}