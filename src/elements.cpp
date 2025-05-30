#include "../include/elements.h"
#include <cmath>     
#include <stdexcept> 

/**
 * @file elements.cpp
 * @brief Implementación de la función para calcular elementos orbitales Keplerianos.
 */


Matrix cross_product_helper( Matrix& v1,  Matrix& v2) {
    if (v1.n_row != 3 || v1.n_column != 1 || v2.n_row != 3 || v2.n_column != 1) {
        throw std::invalid_argument("cross_product_helper: Los vectores deben ser 3x1.");
    }
    Matrix result(3,1);
    result(1,1) = v1(2,1)*v2(3,1) - v1(3,1)*v2(2,1);
    result(2,1) = v1(3,1)*v2(1,1) - v1(1,1)*v2(3,1);
    result(3,1) = v1(1,1)*v2(2,1) - v1(2,1)*v1(1,1); 

    return result;
}

double dot_product_helper( Matrix& v1, Matrix& v2) {
    if (v1.n_row != v2.n_row || v1.n_column != v2.n_column || v1.n_column != 1) {

        throw std::invalid_argument("dot_product_helper: Vectores no compatibles o no son 3x1.");
    }
    double dot_val = 0.0;
    for (int i = 1; i <= v1.n_row; ++i) {
        dot_val += v1(i,1) * v2(i,1);
    }
    return dot_val;
}


static double custom_mod_positive(double value, double divisor) {
    if (divisor == 0.0) return value; // Evitar división por cero
    double result = std::fmod(value, divisor);
    if (result < 0.0) {
        result += divisor;
    }
    return result;
}


OrbitalElements elements( Matrix& y_state) {
    if (y_state.n_row != 6 || y_state.n_column != 1) {
        throw std::invalid_argument("elements: El vector de estado y_state debe ser 6x1.");
    }

    OrbitalElements oe;

    // Extraer posición y velocidad
    Matrix r_vec(3,1), v_vec(3,1);
    r_vec(1,1) = y_state(1,1); r_vec(2,1) = y_state(2,1); r_vec(3,1) = y_state(3,1);
    v_vec(1,1) = y_state(4,1); v_vec(2,1) = y_state(5,1); v_vec(3,1) = y_state(6,1);

    // Vector de momento angular h = r x v
 
    Matrix h_vec = r_vec.cross( v_vec); 

    double mag_h = h_vec.norm();
    oe.p_semilatus_rectum = (mag_h * mag_h) / Const::GM_Earth;
    
    oe.Omega_raan_rad = std::atan2(h_vec(1,1), -h_vec(2,1));
    oe.Omega_raan_rad = custom_mod_positive(oe.Omega_raan_rad, Const::pi2);

  
    double h_xy = std::sqrt(h_vec(1,1) * h_vec(1,1) + h_vec(2,1) * h_vec(2,1));
    oe.i_inclination_rad = std::atan2(h_xy, h_vec(3,1));
   
    double r_cross_h_z_comp = r_vec(1,1)*h_vec(2,1) - r_vec(2,1)*h_vec(1,1);
    oe.u_arg_latitude_rad = std::atan2(r_vec(3,1) * mag_h, -r_cross_h_z_comp ); 
    double R_mag = r_vec.norm();

  
    double v_sq = v_vec.dot(v_vec); 
    double val_2_div_R = 2.0 / R_mag;
    double val_v2_div_GM = v_sq / Const::GM_Earth;
    
    if (std::abs(val_2_div_R - val_v2_div_GM) < 1e-14) { 
             oe.a_semimajor_axis = 1.0 / (val_2_div_R - val_v2_div_GM);
        } else { 
             oe.a_semimajor_axis = 1.0 / (val_2_div_R - val_v2_div_GM); 
        
    } 


    // e*cos(E) y e*sin(E)
    double eCosE = 1.0 - R_mag / oe.a_semimajor_axis; // Válido si a != 0
    double eSinE = r_vec.dot( v_vec) / std::sqrt(Const::GM_Earth * oe.a_semimajor_axis); 
   
    // Excentricidad (e)
    double e_sq = eCosE * eCosE + eSinE * eSinE;
    oe.e_eccentricity = std::sqrt(e_sq);

    // Anomalía Excéntrica (E)
    oe.E_ecc_anomaly_rad = std::atan2(eSinE, eCosE);

    // Anomalía Media (M)
    oe.M_mean_anomaly_rad = oe.E_ecc_anomaly_rad - eSinE; // E - e*sin(E)
    oe.M_mean_anomaly_rad = custom_mod_positive(oe.M_mean_anomaly_rad, Const::pi2);

   
     oe.nu_true_anomaly_rad = std::atan2(std::sqrt(1.0 - e_sq) * eSinE, eCosE - e_sq);
    
    

    // Argumento del Pericentro (omega)
    oe.omega_argp_rad = oe.u_arg_latitude_rad - oe.nu_true_anomaly_rad;
    oe.omega_argp_rad = custom_mod_positive(oe.omega_argp_rad, Const::pi2);

  
    return oe;
}