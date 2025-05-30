#ifndef SIGN_H
#define SIGN_H

/**
 * @file sign.h
 * @brief Declaración de la función para obtener el valor absoluto de 'a' con el signo de 'b'.
 */

/**
 * @brief Devuelve el valor absoluto del primer argumento con el signo del segundo argumento.
 *
 * Si b >= 0.0, devuelve abs(a).
 * Si b < 0.0, devuelve -abs(a).
 *
 * @param a El número del cual se tomará la magnitud (valor absoluto).
 * @param b El número del cual se tomará el signo.
 * @return double El valor absoluto de 'a' con el signo de 'b'.
 */
double sign_(double a, double b);

#endif // SIGN_H