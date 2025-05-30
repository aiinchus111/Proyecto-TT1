#include "..\include\accel.h"

#include "..\include\global.h"      
#include "..\include\SAT_Const.h"            

#include "..\include\iers.h"             
#include "..\include\timediff.h"         
#include "..\include\precc_matrix.h"      
#include "..\include\nut_matrix.h"       
#include "..\include\polematrix.h"      
#include "..\include\gha_matrix.h"       
#include "..\include\mjday_tdb.h"
#include "..\include\jpl_eph_de430.h"    
#include "..\include\accelharmonic.h"   
#include "..\include\AccelPointMass.h"    
#include <stdexcept> 
#include <string>    

/**
 * @file accel.cpp
 * @brief Implementación de la función principal de cálculo de aceleración.
 * @details Esta versión utiliza datos y parámetros globales definidos en el namespace `Global`.
 */

Matrix Accel(double time_from_epoch_sec, 
             Matrix& Y_state) {

    if (Y_state.n_row != 6 || Y_state.n_column != 1) {
        throw std::invalid_argument("Accel: Y_state debe ser un vector columna 6x1.");
    }

    // --- 1. Cálculos de Tiempo ---
  
    double current_Mjd_UTC = AuxParams.Mjd_UTC_epoch + time_from_epoch_sec / 86400.0;

    if (eopdata.n_column == 0) { 
         throw std::runtime_error("Accel: eopdata no inicializada o vacía.");
    }
    // La struct devuelta por IERS debe coincidir con la que usa tu main 
    EopResults eops = IERS(eopdata, current_Mjd_UTC, 'l'); 


    TimeDifferences time_diffs = timediff(eops.UT1_UTC, eops.TAI_UTC);

    double Mjd_UT1 = current_Mjd_UTC + eops.UT1_UTC / 86400.0; 
    double Mjd_TT  = current_Mjd_UTC + time_diffs.TT_UTC / 86400.0; 

    // --- 2. Matrices de Transformación de Coordenadas ---
    Matrix P_prec = PrecMatrix(Const::MJD_J2000, Mjd_TT); 
    Matrix N_nut  = NutMatrix(Mjd_TT);                   
    Matrix T_PN_mat = N_nut * P_prec; 
    
    Matrix GHAm_mat = GHAMatrix(Mjd_UT1);                  
    Matrix PoleM_mat = PoleMatrix(eops.x_pole, eops.y_pole); 
    
    Matrix E_ICRF_to_ITRF = PoleM_mat * GHAm_mat * T_PN_mat; // Matriz E en MATLAB

    // --- 3. Efemérides de Cuerpos Perturbadores ---
    double Mjd_TDB_val = Mjday_TDB(Mjd_TT); 
    
    if (PC.n_column == 0 || PC.n_row == 0) { 
         throw std::runtime_error("Accel: no inicializada o vacía.");
    }

    PlanetPositions planet_pos = JPL_Eph_DE430(Mjd_TDB_val); // JPL_Eph_DE430 usa Global::PC internamente

    // --- 4. Cálculo de Aceleraciones ---
    Matrix r_sat_icrf(3,1);
    r_sat_icrf(1,1) = Y_state(1,1);
    r_sat_icrf(2,1) = Y_state(2,1);
    r_sat_icrf(3,1) = Y_state(3,1);

    
    if (Cnm_coeffs.n_column == 0 || Snm_coeffs.n_column == 0) {
        throw std::runtime_error("Accel: no inicializadas.");
    }
    Matrix a_harmonic = AccelHarmonic(r_sat_icrf, E_ICRF_to_ITRF, 
                                      AuxParams.n_max_gravity, 
                                      AuxParams.m_max_gravity);
                                     
    
    Matrix total_accel = a_harmonic; 

    // Perturbaciones de terceros cuerpos (Sol, Luna)
    if (AuxParams.sun_perturbation) {
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Sun, Const::GM_Sun);
    }
    if (AuxParams.moon_perturbation) {
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Moon, Const::GM_Moon);
    }

    // Perturbaciones planetarias
    if (AuxParams.planets_perturbation) {
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Mercury, Const::GM_Mercury);
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Venus,   Const::GM_Venus);
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Mars,    Const::GM_Mars);
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Jupiter, Const::GM_Jupiter);
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Saturn,  Const::GM_Saturn);
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Uranus,  Const::GM_Uranus);
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Neptune, Const::GM_Neptune);
        total_accel = total_accel + AccelPointMass(r_sat_icrf, planet_pos.r_Pluto,   Const::GM_Pluto);
    }
    
  

    // --- 5. Construir el vector de salida dY ---
    Matrix dY(6,1);
    dY(1,1) = Y_state(4,1); 
    dY(2,1) = Y_state(5,1); 
    dY(3,1) = Y_state(6,1); 
    dY(4,1) = total_accel(1,1); 
    dY(5,1) = total_accel(2,1); 
    dY(6,1) = total_accel(3,1); 
    
    return dY;
}