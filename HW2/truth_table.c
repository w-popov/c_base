/**
 * гр.Д01-134 Попов В.Г.
 * Распечатать таблицы истинности задания №4
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <locale.h>

#define BLUE "\033[1;34m"
#define YELLOW "\033[1;33m"
#define NO "\033[0m"

#define LIMIT_4 4

/**
 * Массив известных результатов.
 * Для импликации: [1, 0, 1, 1],
 * для эквивалентности: [1, 0, 0, 1] и тд.
 */
struct IdentityBits
{
    const uint8_t identity_b[LIMIT_4];
};

/**
 * Структура для хранения двух 
 * булевых значений в виде битов одного байта.
 */
struct BoolValues
{
    uint8_t bita : 1;
    uint8_t bitb : 1;
    uint8_t _ : 6;
};

/**
 * Объединение для заполнения бит BoolValues
 * путем записи числа в byte.
 */
union TruthTable
{
    struct BoolValues values;
    uint8_t byte;
};

/**
 * Тип для функции, которая принимает два (целых)булевых
 * аргументаи возвращает булев результат.
 */
typedef uint8_t(*bool_function)(uint8_t, uint8_t);


/**
 * Импликация: A → B = !A || B
 */
uint8_t implication(uint8_t a, uint8_t b)
{
    return (uint8_t)((!a) || b);
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
uint8_t calculate_expression(uint8_t a, uint8_t b, bool_function func)
{
    return func(a, b);
}

/**
 * Распечатать таблицу истинности (база)
 */
void calculate_identity(struct IdentityBits identity_bits, bool_function b_func)
{
    union TruthTable table;
    for (uint8_t i = 0; i != LIMIT_4; ++i)
    {
        table.byte = i;
        uint8_t a = table.values.bita;
        uint8_t b = table.values.bitb;
        printf("%s| %d | %d |   %d   |%s| %d | %d |           %d          |%s\n",
               BLUE, a, b, identity_bits.identity_b[i],
               YELLOW, a, b, calculate_expression(a, b, b_func), NO);
    }
}

/**
 * Распечатать таблицу истинности для импликации
 */
void print_implication_identity()
{
    struct IdentityBits impl_bits = {.identity_b = {1, 0, 1, 1}};
    printf("\n|------------------- импликация ----------------|\n");
    printf("%s| A | B | A → B |%s| A | B |       (!A || B)      | %s\n", BLUE, YELLOW, NO);
    printf("%s|---|---|-------|%s|---|---|----------------------| %s\n", BLUE, YELLOW, NO);
    calculate_identity(impl_bits, implication);
}

/**
 * Распечатать таблицу истинности для эквивалентности
 */
void print_equivalence_identity()
{
    struct IdentityBits equiv_bits = {.identity_b = {1, 0, 0, 1}};
    printf("\n|--------------- эквивалентность ---------------|\n");
    printf("%s| A | B | A ↔ B |%s| A | B | (A && B)||(!A && !B) | %s\n", BLUE, YELLOW, NO);
    printf("%s|---|---|-------|%s|---|---|----------------------| %s\n", BLUE, YELLOW, NO);
    calculate_identity(equiv_bits, equivalence);
}



int main(void)
{
    setlocale(LC_ALL, "en_US.UTF-8");

    print_implication_identity();
    print_equivalence_identity();

    return EXIT_SUCCESS;
}