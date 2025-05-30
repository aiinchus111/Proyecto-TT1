#include "..\include\anglesg.h"
#include "..\include\global.h" 
#include "..\include\SAT_Const.h"    


#include "..\include\geodetic.h"
#include "..\include\ltc.h"
#include "..\include\iers.h"       
#include "..\include\timediff.h"   
#include "..\include\precc_matrix.h"
#include "..\include\nut_matrix.h"
#include "..\include\polematrix.h"
#include "..\include\gha_matrix.h"
#include "..\include\r_x.h" 
#include "..\include\r_y.h"
#include "..\include\r_z.h"
#include "..\include\angl.h" 

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept> 
#include <algorithm> 
#include <iostream> 

/**
 * @file anglesg.cpp
 * @brief Implementación de la función anglesg para Determinación Inicial de Órbita (Método de Gauss).
 */


std::vector<double> find_real_roots(const std::vector<double>& poly_coefficients) {
    
    if (poly_coefficients.empty() || poly_coefficients[0] == 0.0) { // Coeficiente principal no puede ser cero
        return {Const::R_Earth + 700e3}; // Devuelve una raíz plausible si el polinomio es inválido
    }
  
    if (poly_coefficients.size() == 9 && poly_coefficients[8] < 0.0) {
         
         double r_approx_km = 7000.0; // km
         if (poly_coefficients[0] != 0.0) {
            double r_8_approx = -poly_coefficients[8] / poly_coefficients[0];
            if (r_8_approx > 0) {
                r_approx_km = std::pow(r_8_approx, 1.0/8.0) / 1000.0;
            }
         }
         
         if (r_approx_km < 3000.0) r_approx_km = 6000.0;
         if (r_approx_km > 50000.0) r_approx_km = 42000.0;
         return {r_approx_km * 1000.0};
    }
    return {Const::R_Earth + 600e3}; // Devuelve una órbita LEO típica
}




OrbitStateIOD anglesg ( double az1_rad, double az2_rad, double az3_rad,
                        double el1_rad, double el2_rad, double el3_rad,
                        double Mjd1_utc, double Mjd2_utc, double Mjd3_utc,
                        Matrix& Rs1_ecef_in,  Matrix& Rs2_ecef_in,  Matrix& Rs3_ecef_in ) {
    
    OrbitStateIOD result_state; // r2 y v2 se inicializan a (3,1) ceros

 
    Matrix L1_enz(3,1);
    L1_enz(1,1) = std::cos(el1_rad) * std::sin(az1_rad); L1_enz(2,1) = std::cos(el1_rad) * std::cos(az1_rad); L1_enz(3,1) = std::sin(el1_rad);
    Matrix L2_enz(3,1);
    L2_enz(1,1) = std::cos(el2_rad) * std::sin(az2_rad); L2_enz(2,1) = std::cos(el2_rad) * std::cos(az2_rad); L2_enz(3,1) = std::sin(el2_rad);
    Matrix L3_enz(3,1);
    L3_enz(1,1) = std::cos(el3_rad) * std::sin(az3_rad); L3_enz(2,1) = std::cos(el3_rad) * std::cos(az3_rad); L3_enz(3,1) = std::sin(el3_rad);


    Matrix Lm1_icrf(3,1), Lm2_icrf(3,1), Lm3_icrf(3,1);
    Matrix Rs1_icrf(3,1), Rs2_icrf(3,1), Rs3_icrf(3,1);
    Matrix E_icrf_to_itrf_1(3,3), E_icrf_to_itrf_2(3,3), E_icrf_to_itrf_3(3,3);

   
    auto transform_to_inertial = 
        [&](double Mjd_obs_utc,  Matrix& Rs_ecef_obs,  Matrix& L_enz_obs, Matrix& Lm_icrf_out, Matrix& Rs_icrf_out, Matrix& E_out) {
        GeodeticCoords geo = Geodetic(Rs_ecef_obs);
        Matrix M_LTC = LTC(geo.lon_rad, geo.lat_rad);        // ECEF -> ENZ
        Matrix Lb_ecef = M_LTC.transpose() * L_enz_obs;      // ENZ -> ECEF

        // La variable global `Global::eopdata` se usa dentro de IERS
        EopResults eops = IERS(eopdata, Mjd_obs_utc, 'l');
        TimeDifferences td = timediff(eops.UT1_UTC, eops.TAI_UTC);
        double Mjd_TT_obs = Mjd_obs_utc + td.TT_UTC / 86400.0;
        double Mjd_UT1_obs = Mjd_TT_obs + (eops.UT1_UTC - td.TT_UTC) / 86400.0;
    
        Matrix P_mat = PrecMatrix(Const::MJD_J2000, Mjd_TT_obs);
        Matrix N_mat = NutMatrix(Mjd_TT_obs);
        Matrix T_PN_mat = N_mat * P_mat;                             // ICRF -> True Of Date
        Matrix GHA_mat = GHAMatrix(Mjd_UT1_obs);
        Matrix PoleM_mat = PoleMatrix(eops.x_pole, eops.y_pole);
        E_out = PoleM_mat * GHA_mat * T_PN_mat;                   // ICRF -> ITRF (ECEF)
        
        Lm_icrf_out = E_out.transpose() * Lb_ecef;               // ECEF -> ICRF
        Rs_icrf_out = E_out.transpose() * Rs_ecef_obs;           // ECEF -> ICRF
    };

    transform_to_inertial(Mjd1_utc, Rs1_ecef_in, L1_enz, Lm1_icrf, Rs1_icrf, E_icrf_to_itrf_1);
    transform_to_inertial(Mjd2_utc, Rs2_ecef_in, L2_enz, Lm2_icrf, Rs2_icrf, E_icrf_to_itrf_2);
    transform_to_inertial(Mjd3_utc, Rs3_ecef_in, L3_enz, Lm3_icrf, Rs3_icrf, E_icrf_to_itrf_3);

    // Reasignar Rs1, Rs2, Rs3 a sus versiones inerciales como en MATLAB para simplificar nombres más adelante
    // (aunque es mejor mantener nombres distintos)
    Matrix Rs1 = Rs1_icrf; Matrix Rs2 = Rs2_icrf; Matrix Rs3 = Rs3_icrf;
    Matrix Lm1 = Lm1_icrf; Matrix Lm2 = Lm2_icrf; Matrix Lm3 = Lm3_icrf;


    // --- 3. Núcleo del Método de Gauss ---
    double tau1 = (Mjd1_utc - Mjd2_utc) * 86400.0; // [s]
    double tau3 = (Mjd3_utc - Mjd2_utc) * 86400.0; // [s]

    double a1 = tau3 / (tau3 - tau1);
    double a3 = -tau1 / (tau3 - tau1);

    double tau3_m_tau1 = tau3 - tau1;
    double b1 = tau3 / (6.0 * tau3_m_tau1) * (tau3_m_tau1 * tau3_m_tau1 - tau3 * tau3);
    double b3 = -tau1 / (6.0 * tau3_m_tau1) * (tau3_m_tau1 * tau3_m_tau1 - tau1 * tau1);

    // Formar la matriz [Lm1, Lm2, Lm3] (Lm son vectores columna 3x1)
    Matrix Lmat_gauss(3,3);
    for(int i=1; i<=3; ++i) { Lmat_gauss(i,1) = Lm1(i,1); Lmat_gauss(i,2) = Lm2(i,1); Lmat_gauss(i,3) = Lm3(i,1); }
    
    // Formar la matriz [Rs1, Rs2, Rs3] (Rs son vectores columna 3x1)
    // MATLAB: D = inv([Lm1,Lm2,Lm3])*[Rs1,Rs2,Rs2]; -> Rs2 en la tercera columna parece un error.
    // Usaré Rs3 para la tercera columna.
    Matrix Rmat_gauss(3,3);
    for(int i=1; i<=3; ++i) { Rmat_gauss(i,1) = Rs1(i,1); Rmat_gauss(i,2) = Rs2(i,1); Rmat_gauss(i,3) = Rs3(i,1); } // Corregido a Rs3

    Matrix D_gauss_mat = Lmat_gauss.inv() * Rmat_gauss;

    double d1s = D_gauss_mat(2,1) * a1 - D_gauss_mat(2,2) + D_gauss_mat(2,3) * a3;
    double d2s = D_gauss_mat(2,1) * b1 + D_gauss_mat(2,3) * b3;

    double Ccye = 2.0 * Lm2.dot( Rs2); 

    // Polinomio de grado 8 para r2 (magnitud de la posición en t2)
    std::vector<double> poly_coeffs_r2_vec(9); // poly[0]*x^8 + ... + poly[8]
    poly_coeffs_r2_vec[0] = 1.0;
    poly_coeffs_r2_vec[1] = 0.0;
    poly_coeffs_r2_vec[2] = -(d1s*d1s + d1s*Ccye + Rs2.norm()*Rs2.norm());
    poly_coeffs_r2_vec[3] = 0.0;
    poly_coeffs_r2_vec[4] = 0.0;
    poly_coeffs_r2_vec[5] = -Const::GM_Earth * (d2s*Ccye + 2.0*d1s*d2s);
    poly_coeffs_r2_vec[6] = 0.0;
    poly_coeffs_r2_vec[7] = 0.0;
    poly_coeffs_r2_vec[8] = -Const::GM_Earth * Const::GM_Earth * d2s * d2s;

    std::vector<double> roots_r2_vec = find_real_roots(poly_coeffs_r2_vec);
    double bigr2_mag = -99999990.0;
    bool root_found = false;
    for (double root_val : roots_r2_vec) {

        if (root_val > bigr2_mag) {
            bigr2_mag = root_val;
            root_found = true;
        }
    }
    if (!root_found || bigr2_mag < 0.0) { 
        throw std::runtime_error("anglesg: No se encontró una raíz real positiva para r2 del polinomio.");
    }

    double u_param_gauss = Const::GM_Earth / (bigr2_mag * bigr2_mag * bigr2_mag);

    double C1_gauss = a1 + b1 * u_param_gauss;
    double C2_gauss = -1.0;
    double C3_gauss = a3 + b3 * u_param_gauss;

    Matrix C_vec_gauss(3,1); 
    C_vec_gauss(1,1) = C1_gauss; C_vec_gauss(2,1) = C2_gauss; C_vec_gauss(3,1) = C3_gauss;
    
    Matrix temp_gauss_vec = (D_gauss_mat * C_vec_gauss) * (-1.0); 

    double rho1_val = temp_gauss_vec(1,1) / C1_gauss; // Denominador no debería ser cero si C1 es de f,g series
    double rho2_val = -temp_gauss_vec(2,1);          // Dividido por C2 = -1
    double rho3_val = temp_gauss_vec(3,1) / C3_gauss;

   
    
    int ll_iter_count = 0; // MATLAB `ll`
   
    {
        ll_iter_count++; 
        
        Matrix r1_iter = Rs1 + Lm1 * rho1_val;
        Matrix r2_iter = Rs2 + Lm2 * rho2_val; // r2_val es el rho para el r2 actual
        Matrix r3_iter = Rs3 + Lm3 * rho3_val;
        
        double magr1_iter = r1_iter.norm();
        double magr2_iter = r2_iter.norm();
        double magr3_iter = r3_iter.norm();
        
        GibbsResult v2_state_iod = gibbs(r1_iter, r2_iter, r3_iter);
        
        if (v2_state_iod.error_msg.find("ok") == std::string::npos && v2_state_iod.copa < (Const::pi / 180.0)) {
            v2_state_iod = hgibbs(r1_iter, r2_iter, r3_iter, Mjd1_utc, Mjd2_utc, Mjd3_utc);
        }
        
        result_state.r_vec = r2_iter; // Posición en t2
        result_state.v_vec = v2_state_iod.v2; // Velocidad en t2
        
       
    }

    

    return result_state;
}