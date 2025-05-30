#include "../include/vareqn.h"
#include "../include/global.h"
#include <stdexcept> 

/**
 * @file var_eqn.cpp
 * @brief Implementación de la función VarEqn.
 */

Matrix VarEqn(double time_from_epoch_sec, 
              Matrix& yPhi_augmented_state) {

    if (yPhi_augmented_state.n_row != 42 || yPhi_augmented_state.n_column != 1) {
        throw std::invalid_argument("VarEqn: yPhi_augmented_state debe ser un vector columna 42x1.");
    }

    // --- 1. Cálculos de Tiempo y Transformaciones ---
     if (eopdata.n_column == 0) {
         throw std::runtime_error("VarEqn: eopdata no inicializada o vacía.");
    }
    EopResults eops = IERS(eopdata, AuxParams.Mjd_UTC_epoch, 'l'); // Usa MJD_UTC de la época para EOPs

    TimeDifferences time_diffs = timediff(eops.UT1_UTC, eops.TAI_UTC);

    // Mjd_UT1 para GHAMatrix: basado en Mjd_TT de la época y los deltas de tiempo de esa época
    double Mjd_UT1_for_GHA = AuxParams.Mjd_TT_epoch + (eops.UT1_UTC - time_diffs.TT_UTC) / 86400.0;
                                       
    // Mjd_TT actual para Precesión y Nutación
    double Mjd_TT_current = AuxParams.Mjd_TT_epoch + time_from_epoch_sec / 86400.0;

    Matrix P_prec = PrecMatrix(Const::MJD_J2000, Mjd_TT_current);
    Matrix N_nut  = NutMatrix(Mjd_TT_current);
    Matrix T_PN   = N_nut * P_prec; // ICRF -> True Of Date
    
    Matrix GHAm   = GHAMatrix(Mjd_UT1_for_GHA); 
    Matrix PoleM  = PoleMatrix(eops.x_pole, eops.y_pole);
    
    Matrix E_ICRF_to_ITRF = PoleM * GHAm * T_PN; // ICRF -> ITRF (fijo al cuerpo)

    // --- 2. Extraer Estado y Matriz de Transición de Estados (STM) ---
    Matrix r_vec(3,1);
    Matrix v_vec(3,1);
    for(int i=1; i<=3; ++i) {
        r_vec(i,1) = yPhi_augmented_state(i,1);
        v_vec(i,1) = yPhi_augmented_state(i+3,1);
    }

    Matrix Phi_current(6,6); // STM actual
    for (int j_col = 1; j_col <= 6; ++j_col) { // Itera sobre las columnas de Phi
        for (int i_row = 1; i_row <= 6; ++i_row) { // Itera sobre las filas de Phi
            // yPhi_augmented_state es 1-indexado, y los elementos de Phi empiezan en el índice 7
            Phi_current(i_row, j_col) = yPhi_augmented_state(6 * (j_col - 1) + i_row + 6, 1);
        }
    }

    // --- 3. Calcular Aceleración y Gradiente ---
    // AccelHarmonic y G_AccelHarmonic usan Cnm y Snm de aux_params
    Matrix a_accel = AccelHarmonic(r_vec, E_ICRF_to_ITRF, 
        AuxParams.n_max_gravity, AuxParams.m_max_gravity);

    Matrix G_gradient = G_AccelHarmonic(r_vec, E_ICRF_to_ITRF, 
        AuxParams.n_max_gravity, AuxParams.m_max_gravity);

    // --- 4. Construir la Matriz Jacobiana de la Dinámica A(t) (dfdy) ---
    // A(t) = [ 0_3x3   I_3x3 ]
    //        [ G_3x3   0_3x3 ]
    // donde G = da/dr. Se asume da/dv = 0.
    Matrix dfdy(6,6); // Se inicializa a ceros si el constructor de Matrix lo hace

    // Bloque superior derecho: d(dr/dt)/d(v) = I_3x3
    for (int i=1; i<=3; ++i) dfdy(i, i+3) = 1.0;

    // Bloque inferior izquierdo: d(dv/dt)/d(r) = G_gradient
    for (int i=1; i<=3; ++i) {
        for (int j=1; j<=3; ++j) {
            dfdy(i+3, j) = G_gradient(i,j);
        }
    }
    // Los otros bloques (dv/dr y da/dv) son cero y se asume que Matrix se inicializa a ceros.

    // --- 5. Calcular la Derivada de la Matriz de Transición de Estados ---
    // Phip = dfdy * Phi_current  (dPhi/dt = A(t) * Phi(t))
    Matrix Phip = dfdy * Phi_current;

    // --- 6. Ensamblar el Vector de Salida yPhip (42x1) ---
    Matrix yPhip_out(42,1);

    // Derivadas del vector de estado (primeros 6 elementos)
    // dr/dt = v
    for (int i=1; i<=3; ++i) yPhip_out(i,1) = v_vec(i,1);
    // dv/dt = a
    for (int i=1; i<=3; ++i) yPhip_out(i+3,1) = a_accel(i,1);

    // Derivadas de la STM (siguientes 36 elementos, almacenados por columnas)
    for (int j_col = 1; j_col <= 6; ++j_col) { // Itera sobre las columnas de Phip
        for (int i_row = 1; i_row <= 6; ++i_row) { // Itera sobre las filas de Phip
            yPhip_out(6 * (j_col - 1) + i_row + 6, 1) = Phip(i_row, j_col);
        }
    }
    
    return yPhip_out;
}