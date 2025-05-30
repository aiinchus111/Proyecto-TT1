#include "../include/ltc.h"


/**
 * @file ltc.cpp
 * @brief Implementación de la función LTC.
 */

Matrix LTC(double lon, double lat) {
   
    Matrix Rz_lon = R_z(lon);
    Matrix Ry_neg_lat = R_y(-1.0 * lat);




    Matrix M_intermediate = Ry_neg_lat * Rz_lon;

    
    Matrix M_final(3, 3); 

    for (int j = 1; j <= 3; ++j) { 
        M_final(1, j) = M_intermediate(2, j); // Nueva Fila 1 obtiene elementos de la Antigua Fila 2
        M_final(2, j) = M_intermediate(3, j); // Nueva Fila 2 obtiene elementos de la Antigua Fila 3
        M_final(3, j) = M_intermediate(1, j); // Nueva Fila 3 obtiene elementos de la Antigua Fila 1
    }

    return M_final;
}