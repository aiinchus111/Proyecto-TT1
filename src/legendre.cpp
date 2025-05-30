#include "../include/legendre.h"
#include <cmath>     
#include <stdexcept> 

/**
 * @file legendre_functions.cpp
 * @brief Implementación de la función para calcular polinomios de Legendre y sus derivadas.
 */

LegendreOutput Legendre(int n, int m, double fi) {
    if (m > n) {
        throw std::invalid_argument("Legendre: El orden (m) no puede ser mayor que el grado (n).");
    }

    LegendreOutput result(n, m);
    Matrix& pnm = result.pnm;  
    Matrix& dpnm = result.dpnm;

    double cos_fi = std::cos(fi);
    double sin_fi = std::sin(fi);

  
    pnm(1, 1) = 1.0;    // P_0,0 = 1
    dpnm(1, 1) = 0.0;   // dP_0,0/dfi = 0

    if (m >= 1) { // P_1,1 y dP_1,1 solo se calculan si m >= 1
        pnm(2, 2) = std::sqrt(3.0) * cos_fi;    // P_1,1
        dpnm(2, 2) = -std::sqrt(3.0) * sin_fi;  // dP_1,1/dfi
    } 
   
    for (int i_deg = 2; i_deg <= n; ++i_deg) {
        if (i_deg <= m) { // Solo calcular si el orden actual i_deg está dentro del límite m
            double factor_diag = std::sqrt((2.0 * i_deg + 1.0) / (2.0 * i_deg));
            pnm(i_deg + 1, i_deg + 1) = factor_diag * cos_fi * pnm(i_deg, i_deg);
            dpnm(i_deg + 1, i_deg + 1) = factor_diag * (cos_fi * dpnm(i_deg, i_deg) - sin_fi * pnm(i_deg, i_deg));
        }
    }

   
    for (int i_deg = 1; i_deg <= n; ++i_deg) {
        if (i_deg - 1 <= m && i_deg-1 >=0) { // Solo calcular si el orden actual (i_deg-1) es válido y está dentro del límite m
            double factor_h1 = std::sqrt(2.0 * i_deg + 1.0);
            pnm(i_deg + 1, (i_deg - 1) + 1) = factor_h1 * sin_fi * pnm(i_deg, (i_deg - 1) + 1);
            
             if (i_deg > 0) { 
                 pnm(i_deg + 1, i_deg) = factor_h1 * sin_fi * pnm(i_deg, i_deg);
                 dpnm(i_deg + 1, i_deg) = factor_h1 * (cos_fi * pnm(i_deg, i_deg) + sin_fi * dpnm(i_deg, i_deg));
             }
        }
    }
    
   
    int j_ord = 0;
    int k_loop = 2;
    while (true) {
        if (j_ord > m) { // El orden j_ord no debe exceder m. j_ord+1 es el índice de columna.
            break;
        }
        for (int i_deg = k_loop; i_deg <= n; ++i_deg) {
            if (j_ord > i_deg || j_ord > m ) continue; // orden no puede ser > grado y no debe exceder m

            double common_factor_num = (2.0 * i_deg + 1.0);
            double common_factor_den = static_cast<double>(i_deg - j_ord) * (i_deg + j_ord);
           
            if (common_factor_den == 0) continue; 
            double common_factor = std::sqrt(common_factor_num / common_factor_den);

            double term1_factor_pnm = std::sqrt(2.0 * i_deg - 1.0);
            double term1_pnm = term1_factor_pnm * sin_fi * pnm(i_deg, j_ord + 1);
            
            double term2_factor_num = static_cast<double>(i_deg + j_ord - 1) * (i_deg - j_ord - 1);
            double term2_factor_den = (2.0 * i_deg - 3.0);
            
            double term2_pnm = 0.0;
            if (term2_factor_num >= 0 && term2_factor_den > 0) { // term2_factor_den puede ser 0 si 2i-3 = 0
                term2_pnm = std::sqrt(term2_factor_num / term2_factor_den) * pnm(i_deg - 1, j_ord + 1);
            } else if (term2_factor_num == 0) { // Si el numerador es 0, el término es 0
                 term2_pnm = 0.0;
            }
        

            pnm(i_deg + 1, j_ord + 1) = common_factor * (term1_pnm - term2_pnm);

            // Para dpnm
            double term1_dpnm_part1 = term1_factor_pnm * sin_fi * dpnm(i_deg, j_ord + 1);
            double term1_dpnm_part2 = term1_factor_pnm * cos_fi * pnm(i_deg, j_ord + 1);
            
            double term2_dpnm = 0.0;
             if (term2_factor_num >= 0 && term2_factor_den > 0) {
                term2_dpnm = std::sqrt(term2_factor_num / term2_factor_den) * dpnm(i_deg - 1, j_ord + 1);
            } else if (term2_factor_num == 0) {
                 term2_dpnm = 0.0;
            }

            dpnm(i_deg + 1, j_ord + 1) = common_factor * (term1_dpnm_part1 + term1_dpnm_part2 - term2_dpnm);
        }
        j_ord = j_ord + 1;
        k_loop = k_loop + 1; 
        
        if (j_ord > m) {
             break;
        }
    }
    return result;
}