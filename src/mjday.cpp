#include "../include/mjday.h"
#include <cmath> 

/**
 * @file mjday.cpp
 * @brief Implementación de la función para calcular la Fecha Juliana Modificada (MJD).
 */

double Mjday(int yr, int mon, int day, int hr, int min, double sec) {
    // Algoritmo para calcular el Día Juliano (JD)
    // Los términos .0 aseguran que las operaciones se realicen en punto flotante donde sea necesario.
    double jd = 367.0 * yr
                - std::floor((7.0 * (yr + std::floor((mon + 9.0) / 12.0))) * 0.25)
                + std::floor(275.0 * mon / 9.0)
                + day + 1721013.5;

    // Añadir la fracción del día correspondiente a la hora, minuto y segundo
    jd += ((sec / 60.0 + min) / 60.0 + hr) / 24.0;

    // Convertir Día Juliano (JD) a Fecha Juliana Modificada (MJD)
    double Mjd = jd - 2400000.5;

    return Mjd;
}