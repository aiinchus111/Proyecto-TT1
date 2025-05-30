#include "../include/angl.h"
#include <limits> 

/**
 * @file angl.cpp
 * @brief Implementación de la función para calcular el ángulo entre dos vectores.
 */




double angl(Matrix& vec1, Matrix& vec2) {
    const double small_threshold_sq = 1.0e-16; // (0.00000001)^2
    const double undefined_angle_val = 999999.1;

    double magv1 = vec1.norm(); 
    double magv2 = vec2.norm(); 

    if ((magv1 * magv2) > small_threshold_sq) {

        double dot_v1_v2 = vec1.dot( vec2); 
        
        double temp_cos_theta = dot_v1_v2 / (magv1 * magv2);

        if (temp_cos_theta > 1.0) {
            temp_cos_theta = 1.0;
        } else if (temp_cos_theta < -1.0) {
            temp_cos_theta = -1.0;
        }
       

        return std::acos(temp_cos_theta); // Devuelve ángulo en [0, pi]
    } else {
        return undefined_angle_val;
    }
}