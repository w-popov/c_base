/**
 * Распечатать таблицы истинности задания №4
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <locale.h>



struct BoolValues
{
    uint8_t OO : 1;
    uint8_t OI : 1;
    uint8_t IO : 1;
    uint8_t II : 1;
    uint8_t _ : 4;
};

union TruthTable
{
    struct BoolValues values;
    uint8_t byte;
};

/**
 * Импликация: A → B = !A || B
 */
uint8_t implication(uint8_t a, uint8_t b)
{
    return (uint8_t)(!a || b);
}

/**
 * Эквивалентность: A ↔ B = (A && B) || (!A && !b)
 */
uint8_t equivalence(uint8_t a, uint8_t b)
{
    return (uint8_t)((a && b) || (!a && !b));
}

/**
 * Вычислить выражение, используя заданную функцию
 */
uint8_t calculate_expression(uint8_t a, uint8_t b, uint8_t(*func)(uint8_t, uint8_t))
{
    return func(a, b);
}

int main(void)
{
    setlocale(LC_ALL, "en_US.UTF-8");


    return EXIT_SUCCESS;
}