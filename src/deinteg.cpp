// #include "../include/deinteg.h"
// #include "../include/sign.h" 
// #include <cmath>
// #include <vector>
// #include <numeric>
// #include <algorithm>
// #include <limits>
// #include <functional>


// const double MACHINE_EPSILON_DE = std::numeric_limits<double>::epsilon();
// const double TWOU_DE = 2.0 * MACHINE_EPSILON_DE;
// const double FOURU_DE = 4.0 * MACHINE_EPSILON_DE;

// static const std::vector<double> TWO_POWERS_DE = {
//     1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
//     256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0
// };

// static const std::vector<double> GSTR_COEFFS_DE = {
//     1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,
//     0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
//     0.00592, 0.00524, 0.00468
// };


// Matrix DEInteg(std::function<Matrix(double, const Matrix&)> func,
//                double t, double tout,
//                double relerr, double abserr,
//                int n_eqn, Matrix y_in) {

//     Matrix y = y_in; 
//     DE_Status current_status_internal = DE_Status::DE_INIT; 
//     bool permit_tout_pass = true; 
//     double told_internal = 0.0; 

//     double current_t_internal = t; 
//     Matrix current_yy = y; 

//     Matrix wt(n_eqn, 1);
//     Matrix p_predicted(n_eqn, 1);
//     Matrix yp_derivative(n_eqn, 1);
//     Matrix phi(n_eqn, 17); 

//     std::vector<double> g_coeffs(14);     
//     std::vector<double> sig_coeffs(14);   
//     std::vector<double> rho_coeffs(14);   
//     std::vector<double> w_coeffs(13);     
//     std::vector<double> alpha_coeffs(13); 
//     std::vector<double> beta_coeffs(13);  
//     std::vector<double> v_coeffs(13);     
//     std::vector<double> psi_history(13);  

//     bool start_flag = true; 
//     double h_step = 0.0;    
//     double hold_step = 0.0; 
//     double hnew_step = 0.0; 
//     int k_order = 0;        
//     int kold_order = 0;     
    
//     double delsgn_direction = 1.0; 
//     bool phase1_flag = false;
//     bool nornd_flag = false; 
//     int ifail_counter = 0;      
//     bool crash_flag = false; 
//     bool old_permit_val = true; 


//     if (t == tout) {
//         return y; 
//     }

//     double epsilon_err = std::max(relerr, abserr);

//     if ( (relerr < 0.0) || (abserr < 0.0) || (epsilon_err <= 0.0) ||
//          ( static_cast<int>(current_status_internal) > static_cast<int>(DE_Status::DE_INVPARAM) ) ||
//          ( (current_status_internal != DE_Status::DE_INIT) && (t != told_internal) ) ) {
//         current_status_internal = DE_Status::DE_INVPARAM; 
//         return y; 
//     }

//     double del_t_total = tout - t;
//     double abs_del_t_total = std::abs(del_t_total);

//     double tend_internal = t + 100.0 * del_t_total;
//     if (!permit_tout_pass) {
//         tend_internal = tout;
//     }

//     int num_steps_taken = 0;
//     int kle4_counter = 0;
//     bool stiff_problem = false;
//     double releps_scaled = relerr / epsilon_err;
//     double abseps_scaled = abserr / epsilon_err;

//     if ( (current_status_internal == DE_Status::DE_INIT) || (!old_permit_val) || (delsgn_direction * del_t_total <= 0.0) ) {
//         start_flag = true;
//         current_t_internal = t;
//         current_yy = y;
//         delsgn_direction = sign(1.0, del_t_total);
//         h_step = sign( std::max(FOURU_DE * std::abs(current_t_internal), std::abs(tout - current_t_internal)), tout - current_t_internal );
//     }

//     while (true) { 
//         if (std::abs(current_t_internal - t) >= abs_del_t_total) {
//             Matrix y_out_interpolated(n_eqn, 1);
//             Matrix yp_out_interpolated(n_eqn, 1); 
            
//             g_coeffs[1-1] = 1.0; 
//             rho_coeffs[1-1] = 1.0;
//             double hi_interp = tout - current_t_internal;
//             int ki_interp = kold_order + 1;
            
//             for (int i_idx = 1; i_idx <= ki_interp; ++i_idx) {
//                 w_coeffs[i_idx-1] = 1.0 / static_cast<double>(i_idx);
//             }
            
//             double term_interp = 0.0;
//             for (int j_idx = 2; j_idx <= ki_interp; ++j_idx) {
//                 double psijm1 = psi_history[j_idx-1-1]; 
//                 double gamma_interp = (hi_interp + term_interp) / psijm1;
//                 double eta_interp = hi_interp / psijm1;
//                 for (int i_idx = 1; i_idx <= ki_interp + 1 - j_idx; ++i_idx) {
//                     w_coeffs[i_idx-1] = gamma_interp * w_coeffs[i_idx-1] - eta_interp * w_coeffs[i_idx+1-1];
//                 }
//                 g_coeffs[j_idx-1] = w_coeffs[1-1]; 
//                 rho_coeffs[j_idx-1] = gamma_interp * rho_coeffs[j_idx-1-1];
//                 term_interp = psijm1;
//             }
            
//             for (int j_idx = 1; j_idx <= ki_interp; ++j_idx) {
//                 int i_val = ki_interp + 1 - j_idx;
//                 for(int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                     y_out_interpolated(l_eq,1) = y_out_interpolated(l_eq,1) + g_coeffs[i_val-1] * phi(l_eq, i_val);
//                     yp_out_interpolated(l_eq,1) = yp_out_interpolated(l_eq,1) + rho_coeffs[i_val-1] * phi(l_eq, i_val);
//                 }
//             }
//             y_out_interpolated = current_yy + y_out_interpolated * hi_interp;
//             y = y_out_interpolated;
//             current_status_internal = DE_Status::DE_DONE; 
//             t = tout; 
//             told_internal = t; 
//             old_permit_val = permit_tout_pass;
//             return y; 
//         } 
        
//         if ( !permit_tout_pass && ( std::abs(tout-current_t_internal) < FOURU_DE*std::abs(current_t_internal) ) ) {
//             h_step = tout - current_t_internal;
//             yp_derivative = func(current_t_internal, current_yy); 
//             y = current_yy + yp_derivative * h_step; 
//             current_status_internal = DE_Status::DE_DONE; 
//             t = tout; 
//             told_internal = t; 
//             old_permit_val = permit_tout_pass;
//             return y; 
//         }
        
//         h_step = sign_(std::min(std::abs(h_step), std::abs(tend_internal - current_t_internal)), h_step);
//         for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//             wt(l_eq,1) = releps_scaled * std::abs(current_yy(l_eq,1)) + abseps_scaled;
//         }
        
//         if (std::abs(h_step) < FOURU_DE * std::abs(current_t_internal)) {
//             h_step = sign_(FOURU_DE * std::abs(current_t_internal), h_step);
//             crash_flag = true;
//             current_status_internal = DE_Status::DE_BADACC;
//             relerr = epsilon_err * releps_scaled;
//             abserr = epsilon_err * abseps_scaled;
//             y = current_yy;
//             t = current_t_internal;
//             told_internal = t;
//             old_permit_val = true;
//             return y;
//         }

//         double p5eps = 0.5 * epsilon_err;
//         crash_flag = false;
//         g_coeffs[1-1] = 1.0; 
//         g_coeffs[2-1] = 0.5;
//         sig_coeffs[1-1] = 1.0;

//         ifail_counter = 0;

//         double round_err_est = 0.0;
//         for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//             round_err_est = round_err_est + (y(l_eq,1) * y(l_eq,1)) / (wt(l_eq,1) * wt(l_eq,1));
//         }
//         round_err_est = TWOU_DE * std::sqrt(round_err_est);
//         if (p5eps < round_err_est) {
//             epsilon_err = 2.0 * round_err_est * (1.0 + FOURU_DE);
//             crash_flag = true;
//             current_status_internal = DE_Status::DE_BADACC;
//             relerr = epsilon_err * releps_scaled;
//             abserr = epsilon_err * abseps_scaled;
//             y = current_yy;
//             t = current_t_internal;
//             told_internal = t;
//             old_permit_val = true;
//             return y;
//         }

//         if (start_flag) {
//             yp_derivative = func(current_t_internal, y);
//             double sum_init_deriv = 0.0;
//             for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                 phi(l_eq, 1+1) = yp_derivative(l_eq,1); 
//                 phi(l_eq, 2+1) = 0.0;                   
//                 sum_init_deriv = sum_init_deriv + (yp_derivative(l_eq,1) * yp_derivative(l_eq,1)) / (wt(l_eq,1) * wt(l_eq,1));
//             }
//             sum_init_deriv = std::sqrt(sum_init_deriv);
//             double absh_init = std::abs(h_step);
//             if (epsilon_err < 16.0 * sum_init_deriv * h_step * h_step) {
//                 absh_init = 0.25 * std::sqrt(epsilon_err / sum_init_deriv);
//             }
//             h_step = sign_(std::max(absh_init, FOURU_DE * std::abs(current_t_internal)), h_step);
//             hold_step = 0.0;
//             hnew_step = 0.0;
//             k_order = 1;
//             kold_order = 0;
//             start_flag = false;
//             phase1_flag = true;
//             nornd_flag = true;
//             if (p5eps <= 100.0 * round_err_est) {
//                 nornd_flag = false;
//                 for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                     phi(l_eq, 15+1) = 0.0; 
//                 }
//             }
//         }
        
//         bool step_succeeded_inner;
//         while(true) { 
//             int kp1 = k_order + 1;
//             int kp2 = k_order + 2;
//             int km1 = k_order - 1;
//             int km2 = k_order - 2;
            
//             int ns_steps = 0;
//             if (h_step != hold_step) {
//                 ns_steps = 0;
//             }
//             if (ns_steps <= kold_order) {
//                 ns_steps = ns_steps + 1;
//             }
//             int nsp1 = ns_steps + 1;
            
//             if (k_order >= ns_steps) {
//                 beta_coeffs[ns_steps-1] = 1.0; 
//                 double realns_val = static_cast<double>(ns_steps);
//                 alpha_coeffs[ns_steps-1] = 1.0 / realns_val; 
//                 double temp1_block1 = h_step * realns_val;
//                 sig_coeffs[nsp1-1] = 1.0; 
//                 if (k_order >= nsp1) {
//                     for (int i_idx = nsp1; i_idx <= k_order; ++i_idx) {
//                         int im1_idx = i_idx - 1;
//                         double temp2_block1 = psi_history[im1_idx-1];
//                         psi_history[im1_idx-1] = temp1_block1;
//                         beta_coeffs[i_idx-1] = beta_coeffs[im1_idx-1] * psi_history[im1_idx-1] / temp2_block1;
//                         temp1_block1 = temp2_block1 + h_step;
//                         alpha_coeffs[i_idx-1] = h_step / temp1_block1;
//                         double reali_val = static_cast<double>(i_idx);
//                         sig_coeffs[i_idx+1-1] = reali_val * alpha_coeffs[i_idx-1] * sig_coeffs[i_idx-1];
//                     }
//                 }
//                 psi_history[k_order-1] = temp1_block1;
                
//                 if (ns_steps > 1) {
//                     if (k_order > kold_order) {
//                         double temp4_block1 = static_cast<double>(k_order * kp1);
//                         v_coeffs[k_order-1] = 1.0 / temp4_block1;
//                         int nsm2_val = ns_steps - 2;
//                         for (int j_idx = 1; j_idx <= nsm2_val; ++j_idx) {
//                             int i_val = k_order - j_idx;
//                             v_coeffs[i_val-1] = v_coeffs[i_val-1] - alpha_coeffs[j_idx+1-1] * v_coeffs[i_val+1-1];
//                         }
//                     }
//                     int limit1_block1 = kp1 - ns_steps;
//                     double temp5_block1 = alpha_coeffs[ns_steps-1];
//                     for (int iq_idx = 1; iq_idx <= limit1_block1; ++iq_idx) {
//                         v_coeffs[iq_idx-1] = v_coeffs[iq_idx-1] - temp5_block1 * v_coeffs[iq_idx+1-1];
//                         w_coeffs[iq_idx-1] = v_coeffs[iq_idx-1];
//                     }
//                     g_coeffs[nsp1-1] = w_coeffs[1-1];
//                 } else {
//                     for (int iq_idx = 1; iq_idx <= k_order; ++iq_idx) {
//                         double temp3_block1 = static_cast<double>(iq_idx * (iq_idx + 1));
//                         v_coeffs[iq_idx-1] = 1.0 / temp3_block1;
//                         w_coeffs[iq_idx-1] = v_coeffs[iq_idx-1];
//                     }
//                 }
//                 int nsp2_val = ns_steps + 2;
//                 if (kp1 >= nsp2_val) {
//                     for (int i_idx = nsp2_val; i_idx <= kp1; ++i_idx) {
//                         int limit2_block1 = kp2 - i_idx;
//                         double temp6_block1 = alpha_coeffs[i_idx-1-1]; // alpha(i) in MATLAB
//                         for (int iq_idx = 1; iq_idx <= limit2_block1; ++iq_idx) {
//                             w_coeffs[iq_idx-1] = w_coeffs[iq_idx-1] - temp6_block1 * w_coeffs[iq_idx+1-1];
//                         }
//                         g_coeffs[i_idx-1] = w_coeffs[1-1];
//                     }
//                 }
//             } 
            
//             if (k_order >= nsp1) { 
//                 for (int i_idx = nsp1; i_idx <= k_order; ++i_idx) {
//                     double temp1_block2 = beta_coeffs[i_idx-1];
//                     for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                         phi(l_eq, i_idx) = temp1_block2 * phi(l_eq, i_idx);
//                     }
//                 }
//             }
            
//             for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                 phi(l_eq, kp2) = phi(l_eq, kp1);
//                 phi(l_eq, kp1) = 0.0;
//                 p_predicted(l_eq,1) = 0.0;
//             }
//             for (int j_idx = 1; j_idx <= k_order; ++j_idx) {
//                 int i_val = kp1 - j_idx;
//                 int ip1_val = i_val + 1;
//                 double temp2_block2 = g_coeffs[i_val-1];
//                 for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                     p_predicted(l_eq,1) = p_predicted(l_eq,1) + temp2_block2 * phi(l_eq, i_val);
//                     phi(l_eq, i_val) = phi(l_eq, i_val) + phi(l_eq, ip1_val);
//                 }
//             }
//             if (nornd_flag) {
//                 p_predicted = current_yy + p_predicted * h_step;
//             } else {
//                 for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                     double tau_val = h_step * p_predicted(l_eq,1) - phi(l_eq, 15+1); 
//                     p_predicted(l_eq,1) = current_yy(l_eq,1) + tau_val;
//                     phi(l_eq, 16+1) = (p_predicted(l_eq,1) - current_yy(l_eq,1)) - tau_val;
//                 }
//             }
//             double xold_for_retry = current_t_internal;
//             current_t_internal = current_t_internal + h_step;
//             double absh_current = std::abs(h_step);
//             yp_derivative = func(current_t_internal, p_predicted);
            
//             double erkm2_val = 0.0;
//             double erkm1_val = 0.0;
//             double erk_val = 0.0;
            
//             for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                 double temp3_block2 = 1.0 / wt(l_eq,1);
//                 double temp4_block2 = yp_derivative(l_eq,1) - phi(l_eq, 1+1); 
//                 if (km2 > 0) { 
//                     erkm2_val = erkm2_val + ((phi(l_eq, km1) + temp4_block2) * temp3_block2) *
//                                             ((phi(l_eq, km1) + temp4_block2) * temp3_block2);
//                 }
//                 if (km2 >= 0) { 
//                     erkm1_val = erkm1_val + ((phi(l_eq, k_order) + temp4_block2) * temp3_block2) *
//                                             ((phi(l_eq, k_order) + temp4_block2) * temp3_block2);
//                 }
//                 erk_val = erk_val + (temp4_block2 * temp3_block2) * (temp4_block2 * temp3_block2);
//             }
            
//             if (km2 > 0) { 
//                 erkm2_val = absh_current * sig_coeffs[km1-1] * GSTR_COEFFS_DE[km2-1] * std::sqrt(erkm2_val);
//             }
//             if (km2 >= 0) { 
//                 erkm1_val = absh_current * sig_coeffs[k_order-1] * GSTR_COEFFS_DE[km1-1] * std::sqrt(erkm1_val);
//             }
            
//             double temp5_block2 = absh_current * std::sqrt(erk_val);
//             double err_estimate = temp5_block2 * (g_coeffs[k_order-1] - g_coeffs[kp1-1]);
//             erk_val = temp5_block2 * sig_coeffs[kp1-1] * GSTR_COEFFS_DE[k_order-1];
//             int knew_order = k_order;
            
//             if (km2 > 0) { 
//                 if (std::max(erkm1_val, erkm2_val) <= erk_val) {
//                     knew_order = km1;
//                 }
//             }
//             if (km2 == 0) { 
//                 if (erkm1_val <= 0.5 * erk_val) {
//                     knew_order = km1;
//                 }
//             }
            
//             step_succeeded_inner = (err_estimate <= epsilon_err);
            
//             if (!step_succeeded_inner) {
//                 phase1_flag = false; 
//                 current_t_internal = xold_for_retry;
//                 for (int i_idx = 1; i_idx <= k_order; ++i_idx) {
//                     double temp1_block3 = 1.0 / beta_coeffs[i_idx-1];
//                     int ip1_idx = i_idx + 1;
//                     for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                         phi(l_eq, i_idx) = temp1_block3 * (phi(l_eq, i_idx) - phi(l_eq, ip1_idx));
//                     }
//                 }
//                 if (k_order >= 2) {
//                     for (int i_idx = 2; i_idx <= k_order; ++i_idx) {
//                         psi_history[i_idx-1-1] = psi_history[i_idx-1] - h_step;
//                     }
//                 }
//                 ifail_counter = ifail_counter + 1;
//                 double temp2_block3 = 0.5;
//                 if (ifail_counter > 3) {
//                     if (p5eps < 0.25 * erk_val) {
//                         temp2_block3 = std::sqrt(p5eps / erk_val);
//                     }
//                 }
//                 if (ifail_counter >= 3) {
//                     knew_order = 1;
//                 }
//                 h_step = temp2_block3 * h_step;
//                 k_order = knew_order;
//                 if (std::abs(h_step) < FOURU_DE * std::abs(current_t_internal)) {
//                     crash_flag = true;
//                     h_step = sign_(FOURU_DE * std::abs(current_t_internal), h_step);
//                     epsilon_err = epsilon_err * 2.0;
//                     break; 
//                 }
//             } 
//             if (step_succeeded_inner) {
//                 break; 
//             }
//         } 

//         if (crash_flag) {
//             current_status_internal = DE_Status::DE_BADACC;
//             relerr = epsilon_err * releps_scaled; 
//             abserr = epsilon_err * abseps_scaled; 
//             y = current_yy; 
//             t = current_t_internal;
//             told_internal = t;
//             old_permit_val = true;
//             return y; 
//         }
        
//         kold_order = k_order;
//         hold_step = h_step;
        
//         double temp1_block4 = h_step * g_coeffs[k_order+1-1]; // g(kp1)
//         if (nornd_flag) {
//             for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                 current_yy(l_eq,1) = p_predicted(l_eq,1) + temp1_block4 * (yp_derivative(l_eq,1) - phi(l_eq, 1+1));
//             }
//         } else {
//             for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                 double rho_val = temp1_block4 * (yp_derivative(l_eq,1) - phi(l_eq, 1+1)) - phi(l_eq, 16+1);
//                 current_yy(l_eq,1) = p_predicted(l_eq,1) + rho_val;
//                 phi(l_eq, 15+1) = (current_yy(l_eq,1) - p_predicted(l_eq,1)) - rho_val;
//             }
//         }
//         yp_derivative = func(current_t_internal, current_yy);
        
//         for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//             phi(l_eq, k_order+1) = yp_derivative(l_eq,1) - phi(l_eq, 1+1);
//             phi(l_eq, k_order+2) = phi(l_eq, k_order+1) - phi(l_eq, k_order+2);
//         }
//         for (int i_idx = 1; i_idx <= k_order; ++i_idx) {
//             for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                 phi(l_eq, i_idx) = phi(l_eq, i_idx) + phi(l_eq, k_order+1);
//             }
//         }
        
//         double erkp1_val = 0.0;
//         if ( (knew_order == k_order-1) || (k_order == 12) ) {
//             phase1_flag = false;
//         }
        
//         if (phase1_flag) {
//             k_order = k_order + 1;
           
//         } else {
//             if (knew_order == k_order - 1) {
//                 k_order = k_order - 1;
                
//             } else {
//                  int kp1 = k_order + 1;
//                  if (kp1 <= ns_steps) { 
//                     for (int l_eq = 1; l_eq <= n_eqn; ++l_eq) {
//                         erkp1_val = erkp1_val + (phi(l_eq, kp1+1)/wt(l_eq,1))*(phi(l_eq, kp1+1)/wt(l_eq,1));
//                     }
//                     erkp1_val = std::abs(h_step) * GSTR_COEFFS_DE[kp1-1] * std::sqrt(erkp1_val);
//                     if (k_order > 1) {
//                         if (erkm1_val <= std::min(erk_val, erkp1_val)) {
//                             k_order = k_order - 1; 
                          
//                         } else {
//                             if ( (erkp1_val < erk_val) && (k_order != 12) ) {
//                                 k_order = kp1; 
                              
//                             }
//                         }
//                     } else if (erkp1_val < 0.5 * erk_val) {
//                         k_order = kp1; 
                     
//                     }
//                  } 
//             } 
//         } 
        
//         double erk_for_hnew = 0.0; 
//         if (k_order == kold_order) erk_for_hnew = erk_val;
//         else if (k_order == kold_order -1) erk_for_hnew = erkm1_val;
//         else erk_for_hnew = erkp1_val; // Si el orden aumentó

//         if ( phase1_flag || (p5eps >= erk_for_hnew * TWO_POWERS_DE[k_order+1]) ) { // two(k+2) -> TWO_POWERS_DE[k+1] si 0-idx
//             hnew_step = 2.0 * h_step;
//         } else {
//             if (p5eps < erk_for_hnew) {
//                 double temp2_hnew = static_cast<double>(k_order + 1);
//                 double r_hnew = std::pow(p5eps / erk_for_hnew, (1.0 / temp2_hnew));
//                 hnew_step = std::abs(h_step) * std::max(0.5, std::min(0.9, r_hnew));
//                 hnew_step = sign_(std::max(hnew_step, FOURU_DE * std::abs(current_t_internal)), h_step);
//             } else {
//                 hnew_step = h_step;
//             }
//         }
//         h_step = hnew_step;
        
//         num_steps_taken = num_steps_taken + 1; 
        
//         kle4_counter = kle4_counter + 1;
//         if (kold_order > 4) {
//             kle4_counter = 0;
//         }
//         if (kle4_counter >= 50) {
//             stiff_problem = true;
//         }
//     } 
    
//     return y;
// }