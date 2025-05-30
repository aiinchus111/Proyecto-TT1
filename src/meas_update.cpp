#include "../include/meas_update.h"
#include <stdexcept> 


/**
 * @file meas_update.cpp
 * @brief Implementación de la función MeasUpdate del filtro de Kalman.
 */


MeasUpdateResult MeasUpdate( Matrix& x_current, 
                            Matrix& z_measurement, 
                            Matrix& g_predicted_measurement, 
                            Matrix& s_std_devs, 
                            Matrix& G_jacobian, 
                            Matrix& P_current, 
                            int state_dimension_n) {

    int m = z_measurement.n_row; 


    
    Matrix R_noise_cov(m, m); 
    for (int i = 1; i <= m; ++i) {
        double s_val;
        if (s_std_devs.n_row == m) { // Vector columna
            s_val = s_std_devs(i, 1);
        } else { // Vector fila
            s_val = s_std_devs(1, i);
        }
        R_noise_cov(i, i) = s_val * s_val;
    }

   
    Matrix G_transpose = G_jacobian.transpose();
    Matrix temp_inv = R_noise_cov + G_jacobian * P_current * G_transpose;
    Matrix temp_inv_inverted = temp_inv.inv(); 

    Matrix K_gain = P_current * G_transpose * temp_inv_inverted;

   
    Matrix innovation = z_measurement - g_predicted_measurement;
    Matrix x_updated = x_current + K_gain * innovation;

   
    Matrix I_n = eye(state_dimension_n);
    Matrix P_updated = (I_n - K_gain * G_jacobian) * P_current;
   
    // Empaquetar resultados
    MeasUpdateResult result;
    result.K = K_gain;
    result.x = x_updated;
    result.P = P_updated;
    
    return result;
}