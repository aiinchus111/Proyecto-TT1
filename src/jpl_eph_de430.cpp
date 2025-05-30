#include "../include/jpl_eph_de430.h"
#include "../include/global.h" 
#include "../include/SAT_Const.h"       
#include <cmath>        
#include <vector>        
#include <stdexcept>     
#include <iostream>      

/**
 * @file jpl_eph_de430.cpp
 * @brief Implementación de JPL_Eph_DE430 para calcular posiciones planetarias.
 */



Matrix extract_cheb_coeffs_subset(Matrix& PC_row,
                                  int coeff_block_start_idx, 
                                  int num_coeffs_total_body_coord, 
                                  int num_coeffs_per_sub_interval,
                                  int sub_interval_idx_j) {
    Matrix subset(1, num_coeffs_per_sub_interval); // Vector fila para Cheb3D
    int start_in_block = sub_interval_idx_j * num_coeffs_per_sub_interval + 1; // 1-based
    
    for (int k = 0; k < num_coeffs_per_sub_interval; ++k) {

        subset(1, k + 1) = PC_row(1, coeff_block_start_idx + start_in_block -1 + k);
    }
    return subset;
}


PlanetPositions JPL_Eph_DE430(double Mjd_TDB) {
    PlanetPositions pos_results; 
    if (PC.n_column == 0 || PC.n_row == 0) {
        throw std::runtime_error("JPL_Eph_DE430: PC (datos efemérides JPL) no cargados o vacíos.");
    }

    double JD_TDB = Mjd_TDB + 2400000.5;
    int record_idx = -1;

    for (int i = 1; i <= PC.n_row; ++i) {
        if (JD_TDB >= PC(i, 1) && JD_TDB <= PC(i, 2)) {
            record_idx = i;
            break;
        }
    }

    if (record_idx == -1) {
        throw std::runtime_error("JPL_Eph_DE430: Mjd_TDB (" + std::to_string(Mjd_TDB) + 
                                 ") fuera del rango de las efemerides Global::PC.");
    }


    int num_coeffs_in_record = PC.n_column; 
    Matrix PCtemp_row(1, num_coeffs_in_record);
    for(int c=1; c <= num_coeffs_in_record; ++c) {
        PCtemp_row(1,c) = PC(record_idx, c);
    }
    
    double t1_mjd_interval_start = PCtemp_row(1,1) - 2400000.5; // MJD al inicio del intervalo del registro PCtemp
    double dt_from_interval_start = Mjd_TDB - t1_mjd_interval_start; // Días desde el inicio de t1_mjd_interval_start

    int j_sub_idx; // Índice del sub-intervalo (0-based)
    double Mjd0_sub_interval_start_abs; // Inicio MJD absoluto del sub-intervalo

    std::vector<double> vec_Cx_all, vec_Cy_all, vec_Cz_all;
    Matrix Cx_sub(1,1), Cy_sub(1,1), Cz_sub(1,1); // Se redimensionarán según N_coeffs


    vec_Cx_all.clear(); vec_Cy_all.clear(); vec_Cz_all.clear();
    for(int k=0; k<13; ++k) vec_Cx_all.push_back(PCtemp_row(1, 231+k));
    for(int k=0; k<13; ++k) vec_Cy_all.push_back(PCtemp_row(1, 244+k));
    for(int k=0; k<13; ++k) vec_Cz_all.push_back(PCtemp_row(1, 257+k));
    for(int k=0; k<13; ++k) vec_Cx_all.push_back(PCtemp_row(1, 270+k));
    for(int k=0; k<13; ++k) vec_Cy_all.push_back(PCtemp_row(1, 283+k));
    for(int k=0; k<13; ++k) vec_Cz_all.push_back(PCtemp_row(1, 296+k));

    if (dt_from_interval_start >= 0 && dt_from_interval_start <= 16.0) { j_sub_idx = 0; }
    else if (dt_from_interval_start > 16.0 && dt_from_interval_start <= 32.0) { j_sub_idx = 1; }
    else { throw std::runtime_error("JPL_Eph_DE430: dt fuera de rango para Tierra."); }
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start + 16.0 * j_sub_idx;
    
    Cx_sub = Matrix(1,13); Cy_sub = Matrix(1,13); Cz_sub = Matrix(1,13);
    for(int k=0; k<13; ++k) Cx_sub(1,k+1) = vec_Cx_all[13*j_sub_idx + k];
    for(int k=0; k<13; ++k) Cy_sub(1,k+1) = vec_Cy_all[13*j_sub_idx + k];
    for(int k=0; k<13; ++k) Cz_sub(1,k+1) = vec_Cz_all[13*j_sub_idx + k];
    pos_results.r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 16.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(13 * 8, 0.0); vec_Cy_all.assign(13 * 8, 0.0); vec_Cz_all.assign(13 * 8, 0.0);
    for(int set=0; set<8; ++set) {
        int base_idx = 441 + set * 39;
        for(int k=0; k<13; ++k) vec_Cx_all[set*13+k] = PCtemp_row(1, base_idx+k);
        for(int k=0; k<13; ++k) vec_Cy_all[set*13+k] = PCtemp_row(1, base_idx+13+k);
        for(int k=0; k<13; ++k) vec_Cz_all[set*13+k] = PCtemp_row(1, base_idx+26+k);
    }
    if (dt_from_interval_start >=0 && dt_from_interval_start <=4) {j_sub_idx=0;} 
    else if (dt_from_interval_start > 4 && dt_from_interval_start <=8) {j_sub_idx=1;}
    else if (dt_from_interval_start > 8 && dt_from_interval_start <=12) {j_sub_idx=2;}
    else if (dt_from_interval_start >12 && dt_from_interval_start <=16) {j_sub_idx=3;}
    else if (dt_from_interval_start >16 && dt_from_interval_start <=20) {j_sub_idx=4;}
    else if (dt_from_interval_start >20 && dt_from_interval_start <=24) {j_sub_idx=5;}
    else if (dt_from_interval_start >24 && dt_from_interval_start <=28) {j_sub_idx=6;}
    else if (dt_from_interval_start >28 && dt_from_interval_start <=32) {j_sub_idx=7;}
    else { throw std::runtime_error("JPL_Eph_DE430: dt fuera de rango para Luna."); }
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start + 4.0 * j_sub_idx;
    Cx_sub = Matrix(1,13); Cy_sub = Matrix(1,13); Cz_sub = Matrix(1,13);
    for(int k=0; k<13; ++k) Cx_sub(1,k+1) = vec_Cx_all[13*j_sub_idx + k];
    for(int k=0; k<13; ++k) Cy_sub(1,k+1) = vec_Cy_all[13*j_sub_idx + k];
    for(int k=0; k<13; ++k) Cz_sub(1,k+1) = vec_Cz_all[13*j_sub_idx + k];
    pos_results.r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 4.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(11 * 2, 0.0); vec_Cy_all.assign(11 * 2, 0.0); vec_Cz_all.assign(11 * 2, 0.0);
    for(int set=0; set<2; ++set) {
        int base_idx = 753 + set * 33;
        for(int k=0; k<11; ++k) vec_Cx_all[set*11+k] = PCtemp_row(1, base_idx+k);
        for(int k=0; k<11; ++k) vec_Cy_all[set*11+k] = PCtemp_row(1, base_idx+11+k);
        for(int k=0; k<11; ++k) vec_Cz_all[set*11+k] = PCtemp_row(1, base_idx+22+k);
    }
    if (dt_from_interval_start >= 0 && dt_from_interval_start <= 16.0) { j_sub_idx = 0; }
    else if (dt_from_interval_start > 16.0 && dt_from_interval_start <= 32.0) { j_sub_idx = 1; }
    else { throw std::runtime_error("JPL_Eph_DE430: dt fuera de rango para Sol."); }
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start + 16.0 * j_sub_idx;
    Cx_sub = Matrix(1,11); Cy_sub = Matrix(1,11); Cz_sub = Matrix(1,11);
    for(int k=0; k<11; ++k) Cx_sub(1,k+1) = vec_Cx_all[11*j_sub_idx + k];
    for(int k=0; k<11; ++k) Cy_sub(1,k+1) = vec_Cy_all[11*j_sub_idx + k];
    for(int k=0; k<11; ++k) Cz_sub(1,k+1) = vec_Cz_all[11*j_sub_idx + k];
    pos_results.r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 16.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(14 * 4, 0.0); vec_Cy_all.assign(14 * 4, 0.0); vec_Cz_all.assign(14 * 4, 0.0);
    for(int set=0; set<4; ++set) {
        int base_idx = 3 + set * 42;
        for(int k=0; k<14; ++k) vec_Cx_all[set*14+k] = PCtemp_row(1, base_idx+k);
        for(int k=0; k<14; ++k) vec_Cy_all[set*14+k] = PCtemp_row(1, base_idx+14+k);
        for(int k=0; k<14; ++k) vec_Cz_all[set*14+k] = PCtemp_row(1, base_idx+28+k);
    }
    if (dt_from_interval_start >=0 && dt_from_interval_start <=8) {j_sub_idx=0;}
    else if (dt_from_interval_start >8 && dt_from_interval_start <=16) {j_sub_idx=1;}
    else if (dt_from_interval_start >16 && dt_from_interval_start <=24) {j_sub_idx=2;}
    else if (dt_from_interval_start >24 && dt_from_interval_start <=32) {j_sub_idx=3;}
    else { throw std::runtime_error("JPL_Eph_DE430: dt fuera de rango para Mercurio."); }
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start + 8.0 * j_sub_idx;
    Cx_sub = Matrix(1,14); Cy_sub = Matrix(1,14); Cz_sub = Matrix(1,14);
    for(int k=0; k<14; ++k) Cx_sub(1,k+1) = vec_Cx_all[14*j_sub_idx + k];
    for(int k=0; k<14; ++k) Cy_sub(1,k+1) = vec_Cy_all[14*j_sub_idx + k];
    for(int k=0; k<14; ++k) Cz_sub(1,k+1) = vec_Cz_all[14*j_sub_idx + k];
    pos_results.r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 8.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(10 * 2, 0.0); vec_Cy_all.assign(10 * 2, 0.0); vec_Cz_all.assign(10 * 2, 0.0);
    for(int set=0; set<2; ++set) {
        int base_idx = 171 + set * 30;
        for(int k=0; k<10; ++k) vec_Cx_all[set*10+k] = PCtemp_row(1, base_idx+k);
        for(int k=0; k<10; ++k) vec_Cy_all[set*10+k] = PCtemp_row(1, base_idx+10+k);
        for(int k=0; k<10; ++k) vec_Cz_all[set*10+k] = PCtemp_row(1, base_idx+20+k);
    }
    if (dt_from_interval_start >= 0 && dt_from_interval_start <= 16.0) { j_sub_idx = 0; }
    else if (dt_from_interval_start > 16.0 && dt_from_interval_start <= 32.0) { j_sub_idx = 1; }
    else { throw std::runtime_error("JPL_Eph_DE430: dt fuera de rango para Venus."); }
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start + 16.0 * j_sub_idx;
    Cx_sub = Matrix(1,10); Cy_sub = Matrix(1,10); Cz_sub = Matrix(1,10);
    for(int k=0; k<10; ++k) Cx_sub(1,k+1) = vec_Cx_all[10*j_sub_idx + k];
    for(int k=0; k<10; ++k) Cy_sub(1,k+1) = vec_Cy_all[10*j_sub_idx + k];
    for(int k=0; k<10; ++k) Cz_sub(1,k+1) = vec_Cz_all[10*j_sub_idx + k];
    pos_results.r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 16.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(11, 0.0); vec_Cy_all.assign(11, 0.0); vec_Cz_all.assign(11, 0.0);
    int base_idx_mars = 309;
    for(int k=0; k<11; ++k) vec_Cx_all[k] = PCtemp_row(1, base_idx_mars+k);
    for(int k=0; k<11; ++k) vec_Cy_all[k] = PCtemp_row(1, base_idx_mars+11+k);
    for(int k=0; k<11; ++k) vec_Cz_all[k] = PCtemp_row(1, base_idx_mars+22+k);
    j_sub_idx = 0; 
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start;
    Cx_sub = Matrix(1,11); Cy_sub = Matrix(1,11); Cz_sub = Matrix(1,11);
    for(int k=0; k<11; ++k) {Cx_sub(1,k+1)=vec_Cx_all[k]; Cy_sub(1,k+1)=vec_Cy_all[k]; Cz_sub(1,k+1)=vec_Cz_all[k];}
    pos_results.r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 32.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(8,0.0); vec_Cy_all.assign(8,0.0); vec_Cz_all.assign(8,0.0);
    int base_idx_jupiter = 342;
    for(int k=0; k<8; ++k) vec_Cx_all[k] = PCtemp_row(1, base_idx_jupiter+k);
    for(int k=0; k<8; ++k) vec_Cy_all[k] = PCtemp_row(1, base_idx_jupiter+8+k);
    for(int k=0; k<8; ++k) vec_Cz_all[k] = PCtemp_row(1, base_idx_jupiter+16+k);
    j_sub_idx = 0; Mjd0_sub_interval_start_abs = t1_mjd_interval_start;
    Cx_sub = Matrix(1,8); Cy_sub = Matrix(1,8); Cz_sub = Matrix(1,8);
    for(int k=0; k<8; ++k) {Cx_sub(1,k+1)=vec_Cx_all[k]; Cy_sub(1,k+1)=vec_Cy_all[k]; Cz_sub(1,k+1)=vec_Cz_all[k];}
    pos_results.r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 32.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;
    

    vec_Cx_all.assign(7,0.0); vec_Cy_all.assign(7,0.0); vec_Cz_all.assign(7,0.0);
    int base_idx_saturn = 366;
    for(int k=0; k<7; ++k) vec_Cx_all[k] = PCtemp_row(1, base_idx_saturn+k);
    for(int k=0; k<7; ++k) vec_Cy_all[k] = PCtemp_row(1, base_idx_saturn+7+k);
    for(int k=0; k<7; ++k) vec_Cz_all[k] = PCtemp_row(1, base_idx_saturn+14+k);
    j_sub_idx = 0; Mjd0_sub_interval_start_abs = t1_mjd_interval_start;
    Cx_sub = Matrix(1,7); Cy_sub = Matrix(1,7); Cz_sub = Matrix(1,7);
    for(int k=0; k<7; ++k) {Cx_sub(1,k+1)=vec_Cx_all[k]; Cy_sub(1,k+1)=vec_Cy_all[k]; Cz_sub(1,k+1)=vec_Cz_all[k];}
    pos_results.r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 32.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(6,0.0); vec_Cy_all.assign(6,0.0); vec_Cz_all.assign(6,0.0);
    int base_idx_uranus = 387;
    for(int k=0; k<6; ++k) vec_Cx_all[k] = PCtemp_row(1, base_idx_uranus+k);
    for(int k=0; k<6; ++k) vec_Cy_all[k] = PCtemp_row(1, base_idx_uranus+6+k);
    for(int k=0; k<6; ++k) vec_Cz_all[k] = PCtemp_row(1, base_idx_uranus+12+k);
    j_sub_idx = 0; Mjd0_sub_interval_start_abs = t1_mjd_interval_start;
    Cx_sub = Matrix(1,6); Cy_sub = Matrix(1,6); Cz_sub = Matrix(1,6);
    for(int k=0; k<6; ++k) {Cx_sub(1,k+1)=vec_Cx_all[k]; Cy_sub(1,k+1)=vec_Cy_all[k]; Cz_sub(1,k+1)=vec_Cz_all[k];}
    pos_results.r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 32.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(6,0.0); vec_Cy_all.assign(6,0.0); vec_Cz_all.assign(6,0.0);
    int base_idx_neptune = 405;
    for(int k=0; k<6; ++k) vec_Cx_all[k] = PCtemp_row(1, base_idx_neptune+k);
    for(int k=0; k<6; ++k) vec_Cy_all[k] = PCtemp_row(1, base_idx_neptune+6+k);
    for(int k=0; k<6; ++k) vec_Cz_all[k] = PCtemp_row(1, base_idx_neptune+12+k);
    j_sub_idx = 0; Mjd0_sub_interval_start_abs = t1_mjd_interval_start;

    for(int k=0; k<6; ++k) {Cx_sub(1,k+1)=vec_Cx_all[k]; Cy_sub(1,k+1)=vec_Cy_all[k]; Cz_sub(1,k+1)=vec_Cz_all[k];}
    pos_results.r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 32.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;


    vec_Cx_all.assign(6,0.0); vec_Cy_all.assign(6,0.0); vec_Cz_all.assign(6,0.0);
    int base_idx_pluto = 423;
    for(int k=0; k<6; ++k) vec_Cx_all[k] = PCtemp_row(1, base_idx_pluto+k);
    for(int k=0; k<6; ++k) vec_Cy_all[k] = PCtemp_row(1, base_idx_pluto+6+k);
    for(int k=0; k<6; ++k) vec_Cz_all[k] = PCtemp_row(1, base_idx_pluto+12+k);

    for(int k=0; k<6; ++k) {Cx_sub(1,k+1)=vec_Cx_all[k]; Cy_sub(1,k+1)=vec_Cy_all[k]; Cz_sub(1,k+1)=vec_Cz_all[k];}
    pos_results.r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 32.0, Cx_sub, Cy_sub, Cz_sub) * 1000.0;

  
    vec_Cx_all.assign(10 * 4, 0.0); vec_Cy_all.assign(10 * 4, 0.0); 
    Matrix Cz_Nut_zeros(1,10); // Matriz de ceros para Z
    for(int set=0; set<4; ++set) {
        int base_idx = 819 + set * 20;
        for(int k=0; k<10; ++k) vec_Cx_all[set*10+k] = PCtemp_row(1, base_idx+k);
        for(int k=0; k<10; ++k) vec_Cy_all[set*10+k] = PCtemp_row(1, base_idx+10+k);
    }
    if (dt_from_interval_start >=0 && dt_from_interval_start <=8) {j_sub_idx=0;}
    else if (dt_from_interval_start >8 && dt_from_interval_start <=16) {j_sub_idx=1;}
    else if (dt_from_interval_start >16 && dt_from_interval_start <=24) {j_sub_idx=2;}
    else if (dt_from_interval_start >24 && dt_from_interval_start <=32) {j_sub_idx=3;}
    else { throw std::runtime_error("JPL_Eph_DE430: dt fuera de rango para Nutations."); }
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start + 8.0 * j_sub_idx;
    Cx_sub = Matrix(1,10); Cy_sub = Matrix(1,10);
    for(int k=0; k<10; ++k) Cx_sub(1,k+1) = vec_Cx_all[10*j_sub_idx + k];
    for(int k=0; k<10; ++k) Cy_sub(1,k+1) = vec_Cy_all[10*j_sub_idx + k];
    Matrix Nutations_vec = Cheb3D(Mjd_TDB, 10, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 8.0, Cx_sub, Cy_sub, Cz_Nut_zeros);
   
   
    vec_Cx_all.assign(10 * 4, 0.0); vec_Cy_all.assign(10 * 4, 0.0); vec_Cz_all.assign(10 * 4, 0.0);
    for(int set=0; set<4; ++set) {
        int base_idx = 899 + set * 30;
        for(int k=0; k<10; ++k) vec_Cx_all[set*10+k] = PCtemp_row(1, base_idx+k);
        for(int k=0; k<10; ++k) vec_Cy_all[set*10+k] = PCtemp_row(1, base_idx+10+k);
        for(int k=0; k<10; ++k) vec_Cz_all[set*10+k] = PCtemp_row(1, base_idx+20+k);
    }
   
    if (dt_from_interval_start >=0 && dt_from_interval_start <=8) {j_sub_idx=0;} 
    
    else { throw std::runtime_error("JPL_Eph_DE430: dt fuera de rango para Librations (rechequear j_sub_idx)."); }
    Mjd0_sub_interval_start_abs = t1_mjd_interval_start + 8.0 * j_sub_idx;
  
    Cx_sub = Matrix(1,10); Cy_sub = Matrix(1,10); Cz_sub = Matrix(1,10);
    for(int k=0; k<10; ++k) Cx_sub(1,k+1) = vec_Cx_all[10*j_sub_idx + k];
    for(int k=0; k<10; ++k) Cy_sub(1,k+1) = vec_Cy_all[10*j_sub_idx + k];
    for(int k=0; k<10; ++k) Cz_sub(1,k+1) = vec_Cz_all[10*j_sub_idx + k];
    Matrix Librations_vec = Cheb3D(Mjd_TDB, 10, Mjd0_sub_interval_start_abs, Mjd0_sub_interval_start_abs + 8.0, Cx_sub, Cy_sub, Cz_sub);
    
    
    const double EMRAT1_val = 1.0 / (1.0 + 81.30056907419062); 
    
    pos_results.r_Earth = pos_results.r_Earth - (pos_results.r_Moon * EMRAT1_val);
   
    pos_results.r_Mercury = pos_results.r_Mercury - pos_results.r_Earth;
    pos_results.r_Venus   = pos_results.r_Venus   - pos_results.r_Earth;
    pos_results.r_Mars    = pos_results.r_Mars    - pos_results.r_Earth;
    pos_results.r_Jupiter = pos_results.r_Jupiter - pos_results.r_Earth;
    pos_results.r_Saturn  = pos_results.r_Saturn  - pos_results.r_Earth;
    pos_results.r_Uranus  = pos_results.r_Uranus  - pos_results.r_Earth;
    pos_results.r_Neptune = pos_results.r_Neptune - pos_results.r_Earth;
    pos_results.r_Pluto   = pos_results.r_Pluto   - pos_results.r_Earth;
    pos_results.r_Sun     = pos_results.r_Sun     - pos_results.r_Earth;
    
   

    return pos_results;
}