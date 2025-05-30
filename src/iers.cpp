#include "../include/iers.h"
#include <cmath>      
#include <stdexcept>  
#include <vector>     

/**
 * @file iers.cpp
 * @brief Implementación de la función IERS para obtener parámetros de orientación terrestre.
 */



EopResults IERS(Matrix& eop_data, double Mjd_UTC, char interp_option) {
    EopResults results;

    int target_mjd_floor = static_cast<int>(std::floor(Mjd_UTC));
    int col_idx_i = -1;
    int num_eop_cols = eop_data.n_column; 

    // Encontrar el índice de la columna para el MJD (día entero)
    
    for (int j = 1; j <= num_eop_cols; ++j) {
        if (static_cast<int>(std::floor(eop_data(4, j))) == target_mjd_floor) {
            col_idx_i = j;
            break;
        }
    }

    if (col_idx_i == -1) {
        throw std::runtime_error("IERS: MJD_UTC no encontrado en la tabla eop_data.");
    }

    if (interp_option == 'l') {
        if (col_idx_i + 1 > num_eop_cols) {
            throw std::runtime_error("IERS: No hay datos del día siguiente para interpolación lineal.");
        }

        // Fracción del día para interpolación
        double fixf = Mjd_UTC - std::floor(Mjd_UTC);

        // Interpolar linealmente cada parámetro
        results.x_pole  = eop_data(5, col_idx_i) + (eop_data(5, col_idx_i+1) - eop_data(5, col_idx_i)) * fixf;
        results.y_pole  = eop_data(6, col_idx_i) + (eop_data(6, col_idx_i+1) - eop_data(6, col_idx_i)) * fixf;
        results.UT1_UTC = eop_data(7, col_idx_i) + (eop_data(7, col_idx_i+1) - eop_data(7, col_idx_i)) * fixf;
        results.LOD     = eop_data(8, col_idx_i) + (eop_data(8, col_idx_i+1) - eop_data(8, col_idx_i)) * fixf;
        results.dpsi    = eop_data(9, col_idx_i) + (eop_data(9, col_idx_i+1) - eop_data(9, col_idx_i)) * fixf;
        results.deps    = eop_data(10,col_idx_i) + (eop_data(10,col_idx_i+1)- eop_data(10,col_idx_i)) * fixf;
        results.dx_pole = eop_data(11,col_idx_i) + (eop_data(11,col_idx_i+1)- eop_data(11,col_idx_i)) * fixf;
        results.dy_pole = eop_data(12,col_idx_i) + (eop_data(12,col_idx_i+1)- eop_data(12,col_idx_i)) * fixf;
        
        // TAI_UTC no se interpola, se toma del día actual 
        results.TAI_UTC = eop_data(13, col_idx_i);

        // Conversión a radianes para los parámetros angulares
        results.x_pole  /= Const::Arcs;
        results.y_pole  /= Const::Arcs;
        results.dpsi    /= Const::Arcs;
        results.deps    /= Const::Arcs;
        results.dx_pole /= Const::Arcs;
        results.dy_pole /= Const::Arcs;

    } else {
        results.x_pole  = eop_data(5, col_idx_i) / Const::Arcs;
        results.y_pole  = eop_data(6, col_idx_i) / Const::Arcs;
        results.UT1_UTC = eop_data(7, col_idx_i);
        results.LOD     = eop_data(8, col_idx_i);
        results.dpsi    = eop_data(9, col_idx_i) / Const::Arcs;
        results.deps    = eop_data(10,col_idx_i) / Const::Arcs;
        results.dx_pole = eop_data(11,col_idx_i) / Const::Arcs;
        results.dy_pole = eop_data(12,col_idx_i) / Const::Arcs;
        results.TAI_UTC = eop_data(13,col_idx_i);
    }

    return results;
}