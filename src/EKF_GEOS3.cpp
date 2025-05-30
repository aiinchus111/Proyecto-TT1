
// //
// /** * @file EKF_GEOS3.cpp
//  * @brief Función principal para la Determinación Inicial de Órbita usando Filtro de Kalman Extendido.
//  * @details Este programa implementa un flujo de determinación de órbita para el satélite GEOS-3,
//  * utilizando el método de Gauss para una estimación inicial (actualmente con valores fijos)
//  * y un Filtro de Kalman Extendido para refinar la órbita usando observaciones reales.
//  * Incluye modelos de fuerza como el campo gravitatorio armónico terrestre y perturbaciones
//  * de terceros cuerpos.
//  * @mainpage Programa de Determinación de Órbita EKF_GEOS3

//  */
// //-----------------------------------------------------------------------------

// // --- Standard C++ Headers ---
// #include <iostream> 
// #include <vector>
// #include <string>
// #include <fstream>
// #include <sstream>
// #include <iomanip> 
// #include <cmath>   
// #include <cstdlib> 

// // --- Tus Archivos de Cabecera ---
// #include "../include/matrix.h"    
// #include "../include/SAT_Const.h"     
// #include "../include/global.h" 

// // Funciones de astrodinámica traducidas
// #include "../include/mjday.h"
// #include "../include/position.h"
// #include "../include/ltc.h"
// #include "../include/r_x.h" 
// #include "../include/r_y.h"
// #include "../include/r_z.h"
// #include "../include/gmst.h"
// #include "../include/azelpa.h"     // Devuelve AzElPaData
// #include "../include/timeupdate.h"
// #include "../include/meas_update.h"// Devuelve MeasUpdateData
// #include "../include/vareqn.h"    
// #include "../include/accel.h"      
// #include "../include/deinteg.h"    
// #include "../include/iers.h"       // Devuelve IERSdata
// #include "../include/timediff.h"   // Devuelve timediff_data



// #include "../include/anglesg.h"    
// #include "../include/jpl_eph_de430.h" 


// // --- Funciones Auxiliares para Matrix (si no están en tu clase) ---
// /**
//  * @brief Crea una Matrix 1x1 a partir de un valor escalar.
//  * @param val Valor escalar.
//  * @return Matrix Objeto Matrix de 1x1.
//  */
// Matrix scalar_to_matrix(double val) {
//     Matrix m(1, 1);
//     m(1, 1) = val;
//     return m;
// }

// /**
//  * @brief Crea una Matrix inicializada a ceros.
//  * @param rows Número de filas.
//  * @param cols Número de columnas (igual a rows si no se especifica).
//  * @return Matrix Objeto Matrix con todos los elementos a cero.
//  * @note Asume que el constructor `Matrix(rows, cols)` inicializa a ceros.
//  * Si no, esta función necesitaría un bucle para poner los ceros.
//  */
// Matrix zeros_matrix(int rows, int cols = -1) {
//     if (cols == -1) cols = rows; 
//     Matrix m(rows, cols); 
//     return m;
// }
// /**
//  * @brief Concatena dos vectores fila (matrices 1xN y 1xM) en una matriz 1x(N+M).
//  * @param A Primera matriz fila.
//  * @param B Segunda matriz fila.
//  * @return Matrix Matriz concatenada.
//  * @throws std::runtime_error si A o B no son vectores fila.
//  */
// Matrix concat_row_vectors( Matrix& A,  Matrix& B) {
//     if (A.n_row != 1 || B.n_row != 1) {
//         throw std::runtime_error("concat_row_vectors: Ambas matrices deben ser vectores fila (1xN).");
//     }
//     Matrix C(1, A.n_column + B.n_column);
//     for (int j = 1; j <= A.n_column; ++j) C(1, j) = A(1, j);
//     for (int j = 1; j <= B.n_column; ++j) C(1, A.n_column + j) = B(1, j);
//     return C;
// }


// //------------------------------------------------------------------------------
// // main()
// //------------------------------------------------------------------------------
// /**
// * @brief Función principal del programa de determinación de órbita.
// * @details Realiza los siguientes pasos:
// * 1. Inicializa constantes y datos globales (EOP, coeficientes de gravedad, efemérides JPL).
// * 2. Lee las observaciones del satélite GEOS-3.
// * 3. Realiza una Determinación Inicial de Órbita (IOD) - actualmente con valores fijos.
// * 4. Propaga el estado inicial a una época de referencia del filtro.
// * 5. Ejecuta un Filtro de Kalman Extendido (EKF) iterando sobre las observaciones:
// * a. Propaga el estado y la matriz de transición de estados (STM).
// * b. Calcula la covarianza predicha (Time Update).
// * c. Para cada tipo de medida (Az, El, Dist):
// * i. Predice la medida y calcula el Jacobiano de la medida.
// * ii. Realiza la actualización de la medida (Measurement Update) del estado y la covarianza.
// * 6. Propaga el estado final estimado a la época de la primera observación.
// * 7. Compara con una solución de referencia e imprime los errores.
// * @return 0 si la ejecución es exitosa, 1 en caso de error crítico.
// */
// //------------------------------------------------------------------------------   
// int main(){
//     std::cout << std::fixed << std::setprecision(10); // Formato para salida de doubles

//     // --- 1. Inicialización de Datos Globales ---
//     printf("Inicializando datos globales...\n");
    
//     // --- 2. Leer Observaciones ---
//     printf("Leyendo observaciones de GEOS3.txt...\n");
//     std::vector<ObservationData> observations_vec; // Usar std::vector para leer
//     std::ifstream obs_file_stream("../data/GEOS3.txt"); // Asegúrate que esta ruta es correcta
//     if (!obs_file_stream.is_open()) {
//         std::cerr << "CRITICO: No se pudo abrir GEOS3.txt. Saliendo." << std::endl;
//         return 1;
//     }
//     std::string obs_line_str; // Cambiado de 'line' para evitar conflicto con tu variable global 'line' si existiera
//     double Y1_obs_val, M_obs_val, D_obs_val, h_obs_val, mn_obs_val, s_obs_val, az_obs_val, el_obs_val, Dist_obs_val;
//     int obs_read_count = 0;
//     const int NOBS_MAX_FROM_SCRIPT = 46; // Número de observaciones en el script original

//     while (std::getline(obs_file_stream, obs_line_str) && obs_read_count < NOBS_MAX_FROM_SCRIPT) {
//         if (obs_line_str.length() < 54) { // Chequeo de longitud mínima
//             // std::cerr << "Advertencia: Linea de observacion demasiado corta: " << obs_line_str << std::endl;
//             continue;
//         }
//         try {
//             ObservationData obs_item;
//             // Usar std::stod y std::stoi para parseo más robusto que atof
//             Y1_obs_val = std::stod(obs_line_str.substr(0, 4));
//             M_obs_val  = std::stod(obs_line_str.substr(5, 2));
//             D_obs_val  = std::stod(obs_line_str.substr(8, 2));
//             h_obs_val  = std::stod(obs_line_str.substr(12, 2));
//             mn_obs_val = std::stod(obs_line_str.substr(15, 2));
//             s_obs_val  = std::stod(obs_line_str.substr(18, 6));
//             az_obs_val = std::stod(obs_line_str.substr(25, 8)); // Ajustado a 8 como en MATLAB
//             el_obs_val = std::stod(obs_line_str.substr(35, 7)); // Ajustado a 7 como en MATLAB
//             Dist_obs_val = std::stod(obs_line_str.substr(44, 10));// Ajustado a 10 como en MATLAB

//             obs_item.mjd_utc = Mjday(static_cast<int>(Y1_obs_val), static_cast<int>(M_obs_val), static_cast<int>(D_obs_val), 
//                                      static_cast<int>(h_obs_val), static_cast<int>(mn_obs_val), s_obs_val);
//             obs_item.az_rad  = Const::Rad * az_obs_val;
//             obs_item.el_rad  = Const::Rad * el_obs_val;
//             obs_item.range_m = 1000.0 * Dist_obs_val;
//             observations_vec.push_back(obs_item);
//             obs_read_count++;
//         } catch (const std::exception& e) {
//             std::cerr << "Error parseando observacion: " << obs_line_str << " -> " << e.what() << std::endl;
//         }
//     }
//     obs_file_stream.close();
//     printf("%zu observaciones leidas.\n", observations_vec.size());

//     if (observations_vec.size() < 18) { // El IOD de ejemplo usa la obs 1, 9, 18
//         std::cerr << "CRITICO: No hay suficientes observaciones para el IOD de ejemplo. Se necesitan al menos 18. Saliendo." << std::endl;
//         return 1;
//     }
    
//     Matrix obs_matrix(obs_read_count, 4); 
//     for(size_t i=0; i<observations_vec.size(); ++i){
//         obs_matrix(static_cast<int>(i)+1, 1) = observations_vec[i].mjd_utc;
//         obs_matrix(static_cast<int>(i)+1, 2) = observations_vec[i].az_rad;
//         obs_matrix(static_cast<int>(i)+1, 3) = observations_vec[i].el_rad;
//         obs_matrix(static_cast<int>(i)+1, 4) = observations_vec[i].range_m;
//     }


//     double sigma_range = 92.5;     
//     double sigma_az_rad = 0.0224 * Const::Rad; 
//     double sigma_el_rad = 0.0139 * Const::Rad; 

//     double lat_station_rad = Const::Rad * 21.5748; 
//     double lon_station_rad = Const::Rad * (-158.2706); 
//     double alt_station_m = 300.20;        
    
//     Matrix Rs_station_ecef = Position(lon_station_rad, lat_station_rad, alt_station_m); // Position devuelve 3x1.
   

//     // --- 3. Determinación Inicial de Órbita (IOD) ---
//     double Mjd1_iod = obs_matrix(1,1);  // Usando la Matrix obs_matrix
//     double Mjd2_iod = obs_matrix(9,1);
//     double Mjd3_iod = obs_matrix(18,1);

//     printf("Realizando IOD (usando valores Fijos como en tu EKF_GEOS3.cpp)...\n");
//     Matrix Y0_apr_6x1(6,1);
//     // Usando los valores hardcodeados de tu EKF_GEOS3.cpp
//     double r2_data_fixed[]={6221397.62857869,2867713.77965738,3006155.98509949};
//     double v2_data_fixed[]{4645.04725161807,-2752.21591588205,-7507.99940987033};
//     // Asumiendo que tu Matrix puede ser construida/poblada así:
//     for(int i=0; i<3; ++i) Y0_apr_6x1(i+1,1) = r2_data_fixed[i];
//     for(int i=0; i<3; ++i) Y0_apr_6x1(i+3+1,1) = v2_data_fixed[i];
//     // Si tuvieras `anglesg` implementado:
//     // Y0_apr_6x1 = anglesg(obs_matrix(1,2), obs_matrix(9,2), obs_matrix(18,2),
//     //                      obs_matrix(1,3), obs_matrix(9,3), obs_matrix(18,3),
//     //                      Mjd1_iod, Mjd2_iod, Mjd3_iod,
//     //                      Rs_station_ecef, Rs_station_ecef, Rs_station_ecef);
//     std::cout << "IOD Y0_apr (fijo) [0]: " << Y0_apr_6x1(1,1) << std::endl;


//     // --- 4. Época de Referencia del Filtro y Propagación Inicial ---
//     // Mjd0 en tu C++ main es una fecha específica de 1995, pero los datos GEOS3 son de 1975.
//     // Voy a usar la fecha de la 9na observación (Mjd2_iod) como la época inicial del filtro Mjd0_epoch_filter
//     // y Y0_apr_6x1 como el estado en esa época (sin propagación inicial).
//     double Mjd0_epoch_filter = Mjd2_iod; 
//     printf("Epoca de referencia del filtro Mjd0_epoch_filter: %f\n", Mjd0_epoch_filter);

//     // Configurar Global::AuxParams para la simulación
//     fillAuxParams(Mjd0_epoch_filter, 20, 20, true, true, true); 
//     // Mjd_TT_epoch en Global::AuxParams se actualizará dentro del bucle EKF o antes de llamar a DEInteg.
//     // Calcular Mjd_TT_epoch inicial para Global::AuxParams
//    EopResults eops_epoch_filt = IERS(eopdata, AuxParams.Mjd_UTC_epoch, 'l');
//     TimeDifferences td_epoch_filt = timediff(eops_epoch_filt.UT1_UTC, eops_epoch_filt.TAI_UTC);
//     fillAuxParams(AuxParams.Mjd_UTC_epoch, AuxParams.Mjd_UTC_epoch + td_epoch_filt.TT_UTC / 86400.0);


//     Matrix Y_ekf_current = Y0_apr_6x1; // Estado inicial para el EKF en Mjd0_epoch_filter

//     // --- 5. Configuración del EKF ---
//     Matrix P_cov = zeros_matrix(6,6); 
//     for(int i=1; i<=3; ++i) P_cov(i,i) = 1e8;  
//     for(int i=4; i<=6; ++i) P_cov(i,i) = 1e3;  

//     Matrix LT_station_mat = LTC(lon_station_rad, lat_station_rad); 

//     double t_ekf_relative = 0.0; // Tiempo relativo a Mjd0_epoch_filter

//     printf("Iniciando bucle EKF...\n");
//     // --- 6. Bucle de Mediciones del EKF ---
//     for (size_t i_obs = 0; i_obs < observations_vec.size(); ++i_obs) {
//         const auto& current_obs = observations_vec[i_obs];
//         printf("Procesando observacion %zu en MJD_UTC %f\n", i_obs + 1, current_obs.mjd_utc);

//         double t_old_relative = t_ekf_relative;
//         Matrix Y_old_state = Y_ekf_current;
        
//         double Mjd_UTC_obs_curr = current_obs.mjd_utc;
//         t_ekf_relative = (Mjd_UTC_obs_curr - Mjd0_epoch_filter) * 86400.0; 
//         double dt_prop_ekf_step = t_ekf_relative - t_old_relative;

//         // Actualizar Mjd_UTC y Mjd_TT en AuxParams para el instante actual de la observación
//         EopResults iers_data_iter = IERS(eopdata, Mjd_UTC_obs_curr, 'l');
//         TimeDifferences td_iter = timediff(iers_data_iter.UT1_UTC, iers_data_iter.TAI_UTC);
//         fillAuxParams(Mjd_UTC_obs_curr, Mjd_UTC_obs_curr + td_iter.TT_UTC / 86400.0);
        
//         double Mjd_UT1_obs_curr = AuxParams.Mjd_TT_epoch + (iers_data_iter.UT1_UTC - td_iter.TT_UTC) / 86400.0;

//         // Propagación de Estado y STM
//         Matrix yPhi_augmented_in(42,1);
//         for(int k=1; k<=6; ++k) yPhi_augmented_in(k,1) = Y_old_state(k,1);
//         for(int r_idx=1; r_idx<=6; ++r_idx) { // STM Inicial = Identidad
//             for(int c_idx=1; c_idx<=6; ++c_idx) {
//                 yPhi_augmented_in(6 + (c_idx-1)*6 + r_idx, 1) = (r_idx==c_idx) ? 1.0 : 0.0;
//             }
//         }
        
//         auto varEqnFuncAdapter = [&](double t_sec_var, const Matrix& yPhi_arg_var) -> Matrix {
//             // VarEqn usa AuxParams y Cnm/Snm implícitamente
//             return VarEqn(t_sec_var, yPhi_arg_var); 
//         };
        
//         Matrix yPhi_augmented_out = yPhi_augmented_in; 
//         if (std::abs(dt_prop_ekf_step) > 1e-6) { // Solo propagar si hay tiempo
//              yPhi_augmented_out = DEInteg(varEqnFuncAdapter, 0.0, dt_prop_ekf_step, 1e-13, 1e-6, 42, yPhi_augmented_in);
//         }
        
//         Matrix Phi_propagated(6,6);
//         for (int c = 1; c <= 6; ++c) {
//             for (int r_idx = 1; r_idx <= 6; ++r_idx) {
//                 Phi_propagated(r_idx, c) = yPhi_augmented_out(6 + (c - 1) * 6 + r_idx, 1);
//             }
//         }
//         for(int k=1; k<=6; ++k) Y_ekf_current(k,1) = yPhi_augmented_out(k,1); // Estado propagado

//         // Actualización Temporal de P
//         if (std::abs(dt_prop_ekf_step) > 1e-6) {
//             P_cov = TimeUpdate(P_cov, Phi_propagated); // Usa la sobrecarga sin Qdt explícito
//         }

//         // Predicción de la Medida (Topocéntrica)
//         double theta_gast_current = gmst(Mjd_UT1_obs_curr);
//         Matrix U_earth_rot = R_z(theta_gast_current);
        
//         Matrix r_sat_icrf(3,1);
//         // r_sat_icrf = Y_ekf_current.extract_range(1,3); // Tu C++ usa esto
//         for(int k=1; k<=3; ++k) r_sat_icrf(k,1) = Y_ekf_current(k,1);

//         // s = LT*(U*r - Rs_station_ecef)
//         // r es 3x1, U es 3x3, U*r es 3x1. Rs_station_ecef es 3x1.
//         // U*r - Rs es 3x1. LT_station_mat es 3x3. LT*(...) es 3x1.
//         Matrix s_topo = LT_station_mat * (U_earth_rot * r_sat_icrf - Rs_station_ecef);

//         // --- Actualización de Medida para Azimut ---
//         AzElPaData aep_data1 = AzElPa(s_topo);
//         Matrix dAds_1x3 = aep_data1.dAds;
//         Matrix G_az_jac(1,6); 
//         Matrix dAzdr_1x3 = dAds_1x3 * LT_station_mat * U_earth_rot;
//         for(int k=1; k<=3; ++k) G_az_jac(1,k) = dAzdr_1x3(1,k); // dAz/dv = 0
        
//         MeasUpdateData mu_res1 = MeasUpdate(Y_ekf_current, scalar_to_matrix(current_obs.az_rad), 
//                                             scalar_to_matrix(aep_data1.Az), scalar_to_matrix(sigma_az_rad),
//                                             G_az_jac, P_cov, 6);
//         Y_ekf_current = mu_res1.x; P_cov = mu_res1.P;
        
//         // --- Actualización de Medida para Elevación (recalcular s_topo con Y actualizado) ---
//         for(int k=1; k<=3; ++k) r_sat_icrf(k,1) = Y_ekf_current(k,1);
//         s_topo = LT_station_mat * (U_earth_rot * r_sat_icrf - Rs_station_ecef);
//         AzElPaData aep_data2 = AzElPa(s_topo);
//         Matrix dEds_1x3 = aep_data2.dEds;
//         Matrix G_el_jac(1,6);
//         Matrix dEldr_1x3 = dEds_1x3 * LT_station_mat * U_earth_rot;
//         for(int k=1; k<=3; ++k) G_el_jac(1,k) = dEldr_1x3(1,k);
//         MeasUpdateData mu_res2 = MeasUpdate(Y_ekf_current, scalar_to_matrix(current_obs.el_rad),
//                                             scalar_to_matrix(aep_data2.El), scalar_to_matrix(sigma_el_rad),
//                                             G_el_jac, P_cov, 6);
//         Y_ekf_current = mu_res2.x; P_cov = mu_res2.P;

//         // --- Actualización de Medida para Distancia (recalcular s_topo con Y actualizado) ---
//         for(int k=1; k<=3; ++k) r_sat_icrf(k,1) = Y_ekf_current(k,1);
//         s_topo = LT_station_mat * (U_earth_rot * r_sat_icrf - Rs_station_ecef);
//         double dist_calc = s_topo.norm(); // Asume que s_topo (3x1) tiene un método norm()
//         Matrix dDds_1x3 = zeros_matrix(1,3); // d(Dist)/ds
//         if (dist_calc > 1e-3) {
//             dDds_1x3(1,1) = s_topo(1,1) / dist_calc;
//             dDds_1x3(1,2) = s_topo(2,1) / dist_calc;
//             dDds_1x3(1,3) = s_topo(3,1) / dist_calc;
//         }
//         Matrix G_dist_jac(1,6);
//         Matrix dDdr_1x3 = dDds_1x3 * LT_station_mat * U_earth_rot;
//         for(int k=1; k<=3; ++k) G_dist_jac(1,k) = dDdr_1x3(1,k);
//         MeasUpdateData mu_res3 = MeasUpdate(Y_ekf_current, scalar_to_matrix(current_obs.range_m),
//                                             scalar_to_matrix(dist_calc), scalar_to_matrix(sigma_range),
//                                             G_dist_jac, P_cov, 6);
//         Y_ekf_current = mu_res3.x; P_cov = mu_res3.P;
//     } // Fin del bucle EKF

//     // --- 7. Propagación Final y Resultados ---
//     printf("Propagando estado final EKF a la epoca de la primera observacion...\n");
//     double mjd_utc_ultima_obs = observations_vec.back().mjd_utc;
//     double mjd_utc_primera_obs = observations_vec.front().mjd_utc;

//     // Configurar Global::AuxParams para esta propagación final.
//     // La época de referencia para time_from_epoch_sec será mjd_utc_ultima_obs.
//     fillAuxParams(mjd_utc_ultima_obs, 20,20,true,true,true); // n,m,flags como antes
//     Eopsdifferences =eops_final_prop_main = IERS(eopdata, AuxParams.Mjd_UTC_epoch, 'l');
//     timediff_data td_final_prop_main = timediff(eops_final_prop_main.UT1_UTC, eops_final_prop_main.TAI_UTC);
//     fillAuxParams(AuxParams.Mjd_UTC_epoch, AuxParams.Mjd_UTC_epoch + td_final_prop_main.TT_UTC / 86400.0);


//     double dt_prop_final_ekf_sec = (mjd_utc_primera_obs - mjd_utc_ultima_obs) * 86400.0;
    
//     auto accelFuncFinalAdapter = [&](double t_sec, const Matrix& Y_curr) -> Matrix {
//         return Accel(t_sec, Y_curr); // Accel usa Global::AuxParams
//     };
//     Matrix Y0_estimado_final = DEInteg(accelFuncFinalAdapter, 0.0, dt_prop_final_ekf_sec, 1e-13, 1e-6, 6, Y_ekf_current);
    
//     double vY_true_data[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
//     Matrix Y_true(6,1,vY_true_data); // Asume constructor Matrix(rows,cols,data_array)

//     std::cout << "\n--- Resultados Finales (Estimado vs Verdadero en epoca de primera obs) ---" << std::endl;
//     printf("\nError de Estimacion de Posicion:\n");
//     printf("dX %10.1f [m] (%15.5f)\n",Y0_estimado_final(1,1)-Y_true(1,1), Y0_estimado_final(1,1)-Y_true(1,1));
//     printf("dY %10.1f [m] (%15.5f)\n",Y0_estimado_final(2,1)-Y_true(2,1), Y0_estimado_final(2,1)-Y_true(2,1));
//     printf("dZ %10.1f [m] (%15.5f)\n",Y0_estimado_final(3,1)-Y_true(3,1), Y0_estimado_final(3,1)-Y_true(3,1));
//     printf("\nError de Estimacion de Velocidad:\n");
//     printf("dVx %8.1f [m/s] (%15.5f)\n",Y0_estimado_final(4,1)-Y_true(4,1), Y0_estimado_final(4,1)-Y_true(4,1));
//     printf("dVy %8.1f [m/s] (%15.5f)\n",Y0_estimado_final(5,1)-Y_true(5,1), Y0_estimado_final(5,1)-Y_true(5,1));
//     printf("dVz %8.1f [m/s] (%15.5f)\n",Y0_estimado_final(6,1)-Y_true(6,1), Y0_estimado_final(6,1)-Y_true(6,1));

//     return 0;
// }