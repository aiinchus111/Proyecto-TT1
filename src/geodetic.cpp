#include "../include/geodetic.h"
#include <cmath>     
#include <limits>   
#include <stdexcept> 

/**
 * @file geodetic.cpp
 * @brief Implementación de la función para convertir coordenadas cartesianas a geodésicas.
 */

GeodeticCoords Geodetic( Matrix& r_ecef) {
    if (r_ecef.n_row != 3 || r_ecef.n_column != 1) {
        throw std::invalid_argument("Geodetic: El vector de entrada r_ecef debe ser 3x1.");
    }

    const double R_equ = Const::R_Earth; // Radio ecuatorial
    const double f     = Const::f_Earth;

    const double convergence_criterion = std::numeric_limits<double>::epsilon() * R_equ;
    const double e_sq  = f * (2.0 - f); // Cuadrado de la primera excentricidad (e^2)

    double X = r_ecef(1,1);
    double Y = r_ecef(2,1);
    double Z = r_ecef(3,1);

    double rho_sq = X * X + Y * Y; // Cuadrado de la distancia al eje Z (rho^2)

    // Comprobar validez de los datos de entrada
    if (r_ecef.norm() < convergence_criterion) { 
        throw std::runtime_error("Geodetic: Vector de posicion de entrada invalido (norma cero).");
    }

   
    double dZ = e_sq * Z; // Estimación inicial para dZ
    double dZ_new;
    double SinPhi, N_val, Nh_val;

    int max_iter = 10; // Límite de iteraciones para evitar bucles infinitos
    int iter_count = 0;

    while (iter_count < max_iter) {
        double ZdZ = Z + dZ;
        Nh_val   = std::sqrt(rho_sq + ZdZ * ZdZ); // Nh = N + h
        SinPhi = ZdZ / Nh_val;                   // sin(lat)
        
 
        N_val  = R_equ / std::sqrt(1.0 - e_sq * SinPhi * SinPhi);
        
        dZ_new = N_val * e_sq * SinPhi;

        if (std::abs(dZ - dZ_new) < convergence_criterion) {
            break; // Convergencia alcanzada
        }
        dZ = dZ_new;
        iter_count++;
    }
   

    GeodeticCoords result;
    result.lon_rad = std::atan2(Y, X); // Longitud geodésica
    
    
    result.lat_rad = std::atan2(Z + dZ, std::sqrt(rho_sq));

    result.h_m     = Nh_val - N_val;

    return result;
}