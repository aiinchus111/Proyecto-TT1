#include "../include/azelpa.h"
#include <cmath>  
#include <limits>   

/**
 * @file azelpa.cpp
 * @brief Implementación de la función AzElPa.
 */

AzElPaData AzElPa( Matrix& s) {
    AzElPaData result;

    double s1 = s(1,1); // Componente Este
    double s2 = s(2,1); // Componente Norte
    double s3 = s(3,1); // Componente Zenit

    double rho = std::sqrt(s1 * s1 + s2 * s2);

    // Azimut
    result.Az = std::atan2(s1, s2); // atan2(Este, Norte) -> ángulo desde el Norte, positivo hacia el Este
    if (result.Az < 0.0) {
        result.Az += Const::pi2; // Normalizar a [0, 2*pi)
    }

    // Elevación
    // Manejar el caso rho = 0 para la elevación para evitar división por cero si s3 no es 0.
    if (rho == 0.0) {
        if (s3 > 0.0) {
            result.El = Const::pi / 2.0;
        } else if (s3 < 0.0) {
            result.El = -Const::pi / 2.0;
        } else {
            result.El = 0.0; // Origen, El es 0 por convención o indefinido
        }
    } else {
        result.El = std::atan(s3 / rho);
    }

    // Derivadas Parciales
    // result.dAds y result.dEds ya están dimensionadas 1x3 por el constructor de AzElPaData.

    double rho_sq = rho * rho;
    double dot_s_s = s1 * s1 + s2 * s2 + s3 * s3;

    if (rho == 0.0) { // Si rho es cero, s1 y s2 son cero.
        // dAds: La derivada del azimut no está bien definida en el zenit/nadir.
        // MATLAB daría NaN (0/0). Asignamos NaN explícitamente.
        double nan_val = std::numeric_limits<double>::quiet_NaN();
        result.dAds(1,1) = nan_val; 
        result.dAds(1,2) = nan_val;
        result.dAds(1,3) = 0.0;

        // dEds:
        if (dot_s_s == 0.0) { // Vector s es [0,0,0]
            result.dEds(1,1) = nan_val;
            result.dEds(1,2) = nan_val;
            result.dEds(1,3) = nan_val;
        } else {
            // s1=0, s2=0, rho=0.
            // El = atan(s3/0). El es +/- pi/2.
            // dEl/ds1 = -s1*s3 / (rho*dot_s_s) -> 0/0 -> NaN
            // dEl/ds2 = -s2*s3 / (rho*dot_s_s) -> 0/0 -> NaN
            // dEl/ds3 = rho / dot_s_s -> 0 / s3^2 = 0 (si s3 != 0)
            result.dEds(1,1) = nan_val; 
            result.dEds(1,2) = nan_val; 
            result.dEds(1,3) = 0.0; //  0 / (s3*s3)
        }
    } else {
        result.dAds(1,1) =  s2 / rho_sq;
        result.dAds(1,2) = -s1 / rho_sq;
        result.dAds(1,3) =  0.0;

        if (dot_s_s == 0.0) { // Imposible si rho != 0, a menos que s3 sea NaN/inf, pero bueno proteger.
             double nan_val = std::numeric_limits<double>::quiet_NaN();
             result.dEds(1,1) = nan_val;
             result.dEds(1,2) = nan_val;
             result.dEds(1,3) = nan_val;
        } else {
             result.dEds(1,1) = (-s1 * s3 / rho) / dot_s_s;
             result.dEds(1,2) = (-s2 * s3 / rho) / dot_s_s;
             result.dEds(1,3) =  rho / dot_s_s;
        }
    }

    return result;
}