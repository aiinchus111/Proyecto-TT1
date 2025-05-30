


#include "../include/global.h"

#include "../include/SAT_Const.h"      
#include <iostream>
#include <fstream>   
#include <sstream>   
#include <stdexcept> 

// Param AuxParam;

/**
 * @file global_data.cpp
 * @brief Definición e implementación de variables globales y funciones de inicialización de datos.
 * @details Este archivo contiene las instancias reales de las matrices de datos globales
 * (JPL, gravedad, EOP) y la estructura de parámetros auxiliares. También implementa
 * las funciones para cargar estos datos desde archivos con nombres predefinidos.
 */

    // Definición de las variables globales
    Matrix PC;          ///< Implementación de la matriz de coeficientes JPL DE430.
    Matrix Cnm_coeffs;  ///< Implementación de la matriz de coeficientes Cnm del modelo de gravedad.
    Matrix Snm_coeffs;  ///< Implementación de la matriz de coeficientes Snm del modelo de gravedad.
    Matrix eopdata;     ///< Implementación de la matriz de datos EOP.
    AuxParamGlobal AuxParams; ///< Instancia global de los parámetros auxiliares de la simulación.


    /**
     * @brief Implementación de la carga de coeficientes del modelo de gravedad.
     * @details Carga los coeficientes Cnm y Snm desde "GGM03S.txt" hasta NMAX=180, MMAX=180.
     * @return true si la carga fue exitosa, false en caso contrario.
     */
   
    bool load_gravity_coeffs() {
        const std::string ggm_filename = "../data/GGM03S.txt";
        const int n_max_load = 180; 
        const int m_max_load = 180; 

        std::cout << "INFO: load_gravity_coeffs() cargando " << ggm_filename << std::endl;
        Cnm_coeffs = Matrix(n_max_load + 1, n_max_load + 1); 
        Snm_coeffs = Matrix(n_max_load + 1, n_max_load + 1); 

        std::ifstream file(ggm_filename);
        if (!file.is_open()) {
            std::cerr << "ERROR: No se pudo abrir el archivo de coeficientes de gravedad: " << ggm_filename << std::endl;
            return false;
        }
        std::string line;
        int n_val, m_val;
        double c_val, s_val, sigma_c, sigma_s; // Columnas del archivo GGM03S

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            
            if (line.empty() || line[0] == '%' || line.find_first_not_of(" \t\n\v\f\r") == std::string::npos) continue;

            if (ss >> n_val >> m_val >> c_val >> s_val >> sigma_c >> sigma_s) {
                if (n_val <= n_max_load && m_val <= m_max_load && m_val <= n_val) {
                    
                    if (n_val + 1 <= Cnm_coeffs.n_row && m_val + 1 <= Cnm_coeffs.n_column) {
                         Cnm_coeffs(n_val + 1, m_val + 1) = c_val;
                         Snm_coeffs(n_val + 1, m_val + 1) = s_val;
                    }
                }
            }
        }
        file.close();
        
       
        if (Cnm_coeffs.n_column > 0 && Cnm_coeffs.n_column > 0) {
            if (std::abs(Cnm_coeffs(1,1)) < 1e-9 ) { // Si C00 es cero, ponerlo a 1
                Cnm_coeffs(1,1) = 1.0;
            }
        }
        std::cout << "Coeficientes GGM03S cargados en Cnm_coeffs, Snm_coeffs." << std::endl;
        return true;
    }

    /**
     * @brief Implementación de la carga de datos EOP.
     * @details Carga datos EOP desde "../data/eop19620101.txt" usando `load_eop_data_from_file`.
     * @return true si la carga fue exitosa, false en caso contrario.
     */
  void load_eop_data_from_file() {
    eopdata=zeros(13, 21413);
    FILE* fid = fopen("../data/eop19620101.txt", "r");
   
    for(int i = 1; i <= 21413; i++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &eopdata(1, i), &eopdata(2, i), &eopdata(3, i), &eopdata(4, i),
            &eopdata(5, i), &eopdata(6, i), &eopdata(7, i), &eopdata(8, i),
            &eopdata(9, i), &eopdata(10, i),&eopdata(11, i), &eopdata(12, i),
            &eopdata(13, i));
    }
    fclose(fid);
}
void load_DE430Coeff_data_from_file() {
  
  FILE* fid = fopen("../data/DE430Coeff.txt", "r");
 
  for(int i = 1; i <= 2285; i++) {
    for(int j = 1; j <= 1020; j++) {
      fscanf(fid, "%lf", &PC(i, j));
  }
  fclose(fid);
}
}

    /**
     * @brief Implementación de fillAuxParams (versión detallada).
     * @see global.h
     */
    void fillAuxParams(double mjd_utc, int n, int m, bool sun, bool moon, bool planets) {
        AuxParams.Mjd_UTC_epoch = mjd_utc;
    
        AuxParams.n_max_gravity = n;
        AuxParams.m_max_gravity = m;
        AuxParams.sun_perturbation = sun;
        AuxParams.moon_perturbation = moon;
        AuxParams.planets_perturbation = planets;
        std::cout << "AuxParams actualizados (configuracion detallada).\n";
    }

    /**
     * @brief Implementación de fillAuxParams (actualización de tiempos).
     * @see global.h
     */
    void fillAuxParams(double mjd_utc, double mjd_tt) {
       AuxParams.Mjd_UTC_epoch = mjd_utc;
       AuxParams.Mjd_TT_epoch = mjd_tt;
        std::cout << "AuxParams actualizados (tiempos Mjd_UTC y Mjd_TT).\n";
    }
    Matrix Cnm;
    Matrix Snm;
    
    void GGM03S() {
        
        Cnm = zeros(181, 181);
        Snm = zeros(181, 181);
        
        FILE *fid = fopen("..\\data\\GGM03S.txt", "r");
        if (fid== NULL) {
            printf("Fail open GGM03S.txt file\n");
            exit(EXIT_FAILURE);
        }
    
        double aux;
        for(int n = 0; n <= 180; n++) {
            for(int m = 0; m <= n; m++) {
                fscanf(fid,"%lf %lf %lf %lf %lf %lf",
                    &aux,&aux,&(Cnm(n+1, m+1)),&(Snm(n+1, m+1)),
                    &aux,&aux);
            }
        }
        fclose(fid);
    }    
    void auxparam() {
            AuxParams.Mjd_UTC_epoch = 4.974611635416653e+04;
            AuxParams.n_max_gravity      = 20;
            AuxParams.m_max_gravity     = 20;
            AuxParams.sun_perturbation     = 1;
            AuxParams.moon_perturbation    = 1;
            AuxParams.planets_perturbation = 1;
            AuxParams.Mjd_TT_epoch  = 4.974611706231468e+04;
        }