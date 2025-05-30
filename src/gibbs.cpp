#include "../include/gibbs.h"
#include <cmath>    
#include <stdexcept> 

/**
 * @file gibbs_iod.cpp
 * @brief Implementación del método de Gibbs para determinación de órbita.
 */



/**
 * @brief Calcula el vector unitario de un vector dado.
 * @param vec Vector de entrada (Matrix 3x1).
 * @return Matrix Vector unitario (3x1). Devuelve vector cero si la norma es cero.
 */
Matrix unit_vector_helper( Matrix& vec) {
    if (vec.n_row != 3 || vec.n_column != 1) {
        throw std::invalid_argument("unit_vector_helper: El vector debe ser 3x1.");
    }
    double mag = vec.norm();
    if (std::abs(mag) < 1.0e-12) { // Evitar división por cero
        Matrix zero_vec(3,1); // Devuelve vector cero
        return zero_vec;
    }
    return vec * (1.0 / mag);
}


GibbsResult gibbs( Matrix& r1_in,  Matrix& r2_in,  Matrix& r3_in) {
    GibbsResult result;

    const double small_check = 1.0e-8;

    // Validar dimensiones de entrada
    if (r1_in.n_row!=3 || r1_in.n_column!=1 || 
        r2_in.n_row!=3 || r2_in.n_column!=1 || 
        r3_in.n_row!=3 || r3_in.n_column!=1) {
        result.error_msg = "r1,r2,r3 no 3x1"; 
        return result; 
    }


    double magr1 = r1_in.norm();
    double magr2 = r2_in.norm();
    double magr3 = r3_in.norm();

   
    Matrix p_vec = r2_in.cross( r3_in); // r2 x r3
    Matrix q_vec =r3_in.cross( r1_in); // r3 x r1
    Matrix w_vec = r1_in.cross( r2_in); // r1 x r2

    Matrix pn_unit = unit_vector_helper(p_vec);
    Matrix r1n_unit = unit_vector_helper(r1_in);
    
   
    double dot_pn_r1n = pn_unit.dot(r1n_unit);
    result.copa = std::asin(dot_pn_r1n); // asin devuelve en [-pi/2, pi/2]

    // Límite de coplanaridad (0.01745 rad ~ 1 grado)
    if (std::abs(dot_pn_r1n) > 0.017452406) {
        result.error_msg = "not coplanar "; 
    }
   

    Matrix D_vec = p_vec + q_vec + w_vec;
    double magD = D_vec.norm();

    Matrix N_vec = (p_vec * magr1) + (q_vec * magr2) + (w_vec * magr3); // scalar * vector
    double magN = N_vec.norm();
    
    Matrix Nn_unit = unit_vector_helper(N_vec);
    Matrix Dn_unit = unit_vector_helper(D_vec);

    
    if ((std::abs(magD) < small_check) || (std::abs(magN) < small_check) ||
        (Nn_unit.dot( Dn_unit) < (1.0 - small_check*10) ) ) { 
        if (Nn_unit.dot(Dn_unit) < small_check && !(std::abs(magD) < small_check || std::abs(magN) < small_check)) {
            
        }
        result.error_msg = " impossible  "; 
    } else {
        // Solo calcular si no hay error de "impossible"
        result.theta  = angl(r1_in, r2_in);
        result.theta1 = angl(r2_in, r3_in);

        // Realizar el método de Gibbs para encontrar v2
        double r1mr2 = magr1 - magr2;
        double r3mr1 = magr3 - magr1;
        double r2mr3 = magr2 - magr3;

        Matrix S_vec = (r3_in * r1mr2) + (r2_in * r3mr1) + (r1_in * r2mr3);
        Matrix B_vec = D_vec.cross( r2_in);
        
        double L_param_num = magD * magN;
        if (std::abs(L_param_num) < 1e-14) { // Evitar división por cero
             result.error_msg = "impossible (L0)";
             return result;
        }
        double L_param = std::sqrt(Const::GM_Earth / L_param_num);
        
        if (std::abs(magr2) < 1e-14) { // Evitar división por cero
            result.error_msg = "impossible (r2_0)";
            return result;
        }
        double term_over_r2 = L_param / magr2;
        
        result.v2 = (B_vec * term_over_r2) + (S_vec * L_param);
        // Si no había error previo, se mantiene el "ok"
    }
    
    return result;
}