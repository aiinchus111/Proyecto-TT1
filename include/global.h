
#ifndef GLOBAL_H
#define GLOBAL_H

#include "matrix.h" 
#include <string>
#include <vector>

extern Matrix eopdata;
/**
 * @file global.h
 * @brief Define y declara variables globales y funciones de inicialización para datos de la simulación.
 * @details Este archivo centraliza el acceso a datos como coeficientes de efemérides,
 * modelos de gravedad, parámetros de orientación terrestre (EOP) y parámetros auxiliares
 * de la simulación. La inicialización se gestiona mediante funciones específicas dentro
 * del namespace Global.
 */

/**
 * @struct AuxParamGlobal
 * @brief Estructura para almacenar parámetros auxiliares y de configuración de la simulación.
 * @details Contiene información como las épocas de referencia, grado y orden del modelo
 * de gravedad, y flags para activar/desactivar diversas perturbaciones.
 */
struct AuxParamGlobal {
    double Mjd_UTC_epoch; 
    double Mjd_TT_epoch;  
    
    int n_max_gravity;    
    int m_max_gravity;    
    
    bool sun_perturbation;       
    bool moon_perturbation;     
    bool planets_perturbation;   
    
    /**
     * @brief Constructor por defecto.
     * Inicializa los miembros a valores predeterminados.
     */
    AuxParamGlobal() : Mjd_UTC_epoch(0.0), Mjd_TT_epoch(0.0), 
                       n_max_gravity(0), m_max_gravity(0),
                       sun_perturbation(false), moon_perturbation(false), 
                       planets_perturbation(false) {}
};

    
    extern Matrix PC;         
    extern Matrix Cnm_coeffs; 
    extern Matrix Snm_coeffs; 
    extern Matrix eopdata;     
    extern AuxParamGlobal AuxParams; 

    /**
     * @brief Carga los coeficientes de las efemérides JPL DE430.
     * @details Lee los coeficientes desde un archivo predefinido (ej. "../data/DE430Coeff.dat")
     * y los almacena en la variable global `Global::PC`.
     */
    void load_DE430Coeff_data_from_file();

    /**
     * @brief Carga los coeficientes del modelo de gravedad (ej. GGM03S).
     * @details Lee los coeficientes $C_{nm}$ y $S_{nm}$ desde un archivo predefinido
     * (ej. "../data/GGM03S.txt") hasta un grado y orden máximos también predefinidos
     * internamente, y los almacena en `Global::Cnm_coeffs` y `Global::Snm_coeffs`.
     * @return true si la carga fue exitosa, false en caso contrario.
     */
    bool load_gravity_coeffs();

    /**
     * @brief Carga los datos de Parámetros de Orientación Terrestre (EOP).
     * @details Lee los datos EOP desde un archivo predefinido (ej. "../data/eop19620101.txt")
     * y los almacena en la variable global `Global::eopdata`. Utiliza la función
     * `load_eop_data_from_file` (de `eop_loader.h`).
     */
    void load_eop_data_from_file();

    /**
     * @brief Rellena los campos principales de la estructura global `Global::AuxParams`.
     * @param mjd_utc MJD en UTC para la época de referencia.
     * @param n Grado máximo del modelo de gravedad.
     * @param m Orden máximo del modelo de gravedad.
     * @param sun Flag para activar perturbación solar (true=activado).
     * @param moon Flag para activar perturbación lunar (true=activado).
     * @param planets Flag para activar perturbaciones planetarias (true=activado).
     * @note Mjd_TT_epoch no se establece aquí; se actualiza por separado o con la otra sobrecarga.
     */
    void fillAuxParams(double mjd_utc, int n, int m, bool sun, bool moon, bool planets);

    /**
     * @brief Actualiza los tiempos de época (MJD_UTC y MJD_TT) en `Global::AuxParams`.
     * @param mjd_utc MJD en UTC para la nueva época de referencia.
     * @param mjd_tt MJD en TT para la nueva época de referencia.
     */
    void fillAuxParams(double mjd_utc, double mjd_tt); 

    void GGM03S();
    void auxparam();
#endif // GLOBAL_H