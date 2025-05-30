#ifndef MJDAY_H
#define MJDAY_H

/**
 * @file mjday.h
 * @brief Declaración de la función para calcular la Fecha Juliana Modificada (MJD).
 */

/**
 * @brief Calcula la Fecha Juliana Modificada (MJD) a partir de una fecha y hora dadas.
 *
 * La función toma el año, mes, día, hora, minuto y segundo y devuelve
 * la Fecha Juliana Modificada correspondiente.
 * Si la hora, minuto y segundo no se proporcionan, se asumen como 00:00:00.
 *
 * @param yr Año (ej. 2000).
 * @param mon Mes (1-12).
 * @param day Día del mes (1-31).
 * @param hr Hora en tiempo universal (0-23). Valor por defecto es 0.
 * @param min Minuto en tiempo universal (0-59). Valor por defecto es 0.
 * @param sec Segundo en tiempo universal (0-59.999...). Valor por defecto es 0.0.
 * @return double La Fecha Juliana Modificada.
 *
 * @note La fórmula para el cálculo del Día Juliano está basada en algoritmos astronómicos estándar.
 */
double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);

#endif // MJDAY_H