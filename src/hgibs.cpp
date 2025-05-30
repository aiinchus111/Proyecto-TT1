#include "../include/hgibs.h" 
#include <cmath>        
#include <stdexcept>    
#include <vector>      


/**
 * @file hgibbs_iod.cpp
 * @brief Implementación del método de Herrick-Gibbs para determinación de órbita.
 */

GibbsResult hgibbs (Matrix& r1_in,  Matrix& r2_in, Matrix& r3_in,
                    double Mjd1, double Mjd2, double Mjd3) {
    
    GibbsResult result; // error_msg se inicializa a "ok" por defecto

    // Validar dimensiones de entrada
    if (r1_in.n_row!=3 || r1_in.n_column!=1 || 
        r2_in.n_row!=3 || r2_in.n_column!=1 || 
        r3_in.n_row!=3 || r3_in.n_column!=1) {
        result.error_msg = "r1,r2,r3 no 3x1";
        return result; 
    }
    
    // Constantes

    const double GM_Earth_val = Const::GM_Earth; // Usar de const.h
    const double tolangle_rad = 0.01745329251994; // 1 grado en radianes

    // Magnitudes
    double magr1 = r1_in.norm();
    double magr2 = r2_in.norm();
    double magr3 = r3_in.norm();

   
    if (magr1 < 1e-6 || magr2 < 1e-6 || magr3 < 1e-6) {
        result.error_msg = "r magnitud cero";
        return result;
    }

   
    // Intervalos de tiempo en segundos
    double dt21 = (Mjd2 - Mjd1) * 86400.0;
    double dt31 = (Mjd3 - Mjd1) * 86400.0;
    double dt32 = (Mjd3 - Mjd2) * 86400.0;

    // Comprobar si los tiempos son distintos para evitar división por cero
    if (std::abs(dt21) < 1e-6 || std::abs(dt31) < 1e-6 || std::abs(dt32) < 1e-6) {
        result.error_msg = "tiempos iguales";
        return result;
    }
    
    // Comprobación de coplanaridad
    Matrix p_cross_r2r3 = r2_in.cross( r3_in); 
    Matrix pn_unit = unit_vector_helper(p_cross_r2r3);
    Matrix r1n_unit = unit_vector_helper(r1_in);
    
    double dot_pn_r1n_val = pn_unit.dot( r1n_unit); 
    result.copa = std::asin(dot_pn_r1n_val);

    if (std::abs(dot_pn_r1n_val) > 0.017452406) { 
        result.error_msg = "not coplanar "; 
    }

    result.theta  = angl(r1_in, r2_in);
    result.theta1 = angl(r2_in, r3_in);

    if ((result.theta > tolangle_rad) || (result.theta1 > tolangle_rad)) {
   
        result.error_msg = "   angl > 1deg";
    }


    double term_dt21_dt31 = dt21 * dt31;
    double term_dt21_dt32 = dt21 * dt32;
    double term_dt32_dt31 = dt32 * dt31;

    if (std::abs(term_dt21_dt31) < 1e-9 || std::abs(term_dt21_dt32) < 1e-9 || std::abs(term_dt32_dt31) < 1e-9) {
        result.error_msg = "dt prod cero  "; 
        return result;
    }
    
    double term1_coeff = -dt32 * (1.0 / term_dt21_dt31 + GM_Earth_val / (12.0 * magr1 * magr1 * magr1));
    double term2_coeff = (dt32 - dt21) * (1.0 / term_dt21_dt32 + GM_Earth_val / (12.0 * magr2 * magr2 * magr2));
    double term3_coeff = dt21 * (1.0 / term_dt32_dt31 + GM_Earth_val / (12.0 * magr3 * magr3 * magr3));

   
    result.v2 = (r1_in * term1_coeff) + (r2_in * term2_coeff) + (r3_in * term3_coeff);
    
    
    return result;
}