#include "../include/position.h"

#include <cmath>   

/**
 * @file position.cpp
 * @brief Implementación de la función para convertir coordenadas geodésicas a cartesianas.
 */

Matrix Position(double lon, double lat, double h) {
    double R_equ = Const::R_Earth; // Radio ecuatorial de la Tierra [m]
    double f     = Const::f_Earth; // Aplanamiento de la Tierra

    double e2     = f * (2.0 - f);             // Cuadrado de la primera excentricidad
    double CosLat = std::cos(lat);             // Coseno de la latitud geodésica
    double SinLat = std::sin(lat);             // Seno de la latitud geodésica

    // Radio de curvatura en la vertical primera (N)
    double N = R_equ / std::sqrt(1.0 - e2 * SinLat * SinLat);

    
    Matrix r(3, 1); // Vector columna 3x1

    // x = (N + h) * cos(lat) * cos(lon)
    r(1, 1) = (N + h) * CosLat * std::cos(lon);
    // y = (N + h) * cos(lat) * sin(lon)
    r(2, 1) = (N + h) * CosLat * std::sin(lon);
    // z = ((1 - e^2) * N + h) * sin(lat)
    r(3, 1) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}