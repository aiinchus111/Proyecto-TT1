#include "../include/accelharmonic.h"
#include "../include/SAT_Const.h"  

#include <cmath>    
#include <stdexcept>

/**
 * @file accel_harmonic.cpp
 * @brief Implementación de la función AccelHarmonic.
 */



Matrix AccelHarmonic( Matrix& r_inertial, Matrix& E,
                     int n_max, int m_max) {

    const double r_ref = 6378.1363e3;  // Radio de referencia [m]; GGM03S
    const double gm    = 398600.4415e9; // Parámetro gravitacional [m^3/s^2]; GGM03S

    if (m_max > n_max) {
        throw std::invalid_argument("AccelHarmonic: m_max no puede ser mayor que n_max.");
    }
    if (r_inertial.n_row != 3 || r_inertial.n_column != 1 ||
        E.n_row != 3 || E.n_column != 3) {
        r_inertial=r_inertial.transpose();
        E=E.transpose();
    }
    


    Matrix r_bf = E * r_inertial;


    double d = r_bf.norm(); 

    if (d < 1.0e-6) { // Si está prácticamente en el centro
        Matrix a_zero(3,1); 
        return a_zero;
    }

    double r_bf_x = r_bf(1,1);
    double r_bf_y = r_bf(2,1);
    
    double r_bf_z = r_bf(3,1);

    double latgc = std::asin(r_bf_z / d);
    double lon   = std::atan2(r_bf_y, r_bf_x);
    
    LegendreOutput legendre_data = Legendre(n_max, m_max, latgc);
    Matrix& pnm = legendre_data.pnm;
    Matrix& dpnm = legendre_data.dpnm;

    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;
    
    double r_ref_div_d = r_ref / d;
    double r_ref_div_d_pow_n = 1.0; // Para n=0 -> (r_ref/d)^0 = 1

    for (int n = 0; n <= n_max; ++n) {
        if (n > 0) {
            r_ref_div_d_pow_n *= r_ref_div_d; // Actualiza (r_ref/d)^n
        }

        double b1 = (-gm / (d * d)) * r_ref_div_d_pow_n * (n + 1.0);
        double b2 = (gm / d) * r_ref_div_d_pow_n;
        

        double q1_sum_m = 0.0; 
        double q2_sum_m = 0.0; 
        double q3_sum_m = 0.0; 

        for (int m = 0; m <= m_max; ++m) {
            if (m > n) continue; 

            double cos_m_lon = std::cos(static_cast<double>(m) * lon);
            double sin_m_lon = std::sin(static_cast<double>(m) * lon);
            cout<<m<<","<<n;
            double Cnm_val = Cnm_coeffs(n + 1, m + 1);
            double Snm_val = Snm_coeffs(n + 1, m + 1);
            double pnm_val = pnm(n + 1, m + 1);
            double dpnm_val = dpnm(n + 1, m + 1);
            
            q1_sum_m += pnm_val  * (Cnm_val * cos_m_lon + Snm_val * sin_m_lon);
            q2_sum_m += dpnm_val * (Cnm_val * cos_m_lon + Snm_val * sin_m_lon);
            if (m > 0) { 
                 q3_sum_m += static_cast<double>(m) * pnm_val * (Snm_val * cos_m_lon - Cnm_val * sin_m_lon);
            }
        }
        dUdr     += q1_sum_m * b1;
        dUdlatgc += q2_sum_m * b2; 
        dUdlon   += q3_sum_m * b2; // usando b2 para b3
    }

    Matrix a_bf(3,1);
    cout<<149;
    double r2xy = r_bf_x * r_bf_x + r_bf_y * r_bf_y;
    double sqrt_r2xy = std::sqrt(r2xy);
    cout<<149;
    if (sqrt_r2xy < 1.0e-10) { // Cerca del polo
        a_bf(1,1) = (1.0/d * dUdr) * r_bf_x; 
        a_bf(2,1) = (1.0/d * dUdr) * r_bf_y; 
        a_bf(3,1) = (1.0/d * dUdr) * r_bf_z + (sqrt_r2xy / (d*d)) * dUdlatgc; 
    } else {
        double term_dUdr_div_d = dUdr / d;
        double term_dlat_common = dUdlatgc / (d * d * sqrt_r2xy);
        double term_dlon_common = dUdlon / r2xy;

        a_bf(1,1) = (term_dUdr_div_d - r_bf_z * term_dlat_common) * r_bf_x - term_dlon_common * r_bf_y;
        a_bf(2,1) = (term_dUdr_div_d - r_bf_z * term_dlat_common) * r_bf_y + term_dlon_common * r_bf_x;
        a_bf(3,1) = term_dUdr_div_d * r_bf_z + sqrt_r2xy / (d*d) * dUdlatgc;
    }

    Matrix E_transpose = E.transpose(); 
    Matrix a_inertial = E_transpose * a_bf;

    return a_inertial;
}