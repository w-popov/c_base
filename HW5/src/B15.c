/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B15: Дана последовательность ненулевых целых чисел, в конце
 *      последовательности число 0. Посчитать количество четных чисел
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "HW5.h"

typedef unsigned uint;

#ifndef TEST_DEF_HW5
int main (void)
{
    enum { SIZE_ARRAY = INT16_MAX };
    char input[SIZE_ARRAY] = {'\0'};

    // Читать все символы пока не встретится \n
    scanf(" %[^\n]", input);
    printf("%u\n", number_of_even_numbers(input, SIZE_ARRAY));

    return EXIT_SUCCESS;
}
#endif

/**
 * Вернуть количество четных чисел
 */
uint number_of_even_numbers (char *input_numbers, const int size)
{
    enum { DIVIDER = 2 };
    uint counter_even = 0;
    uint number = UINT32_MAX, offset = 0;
    int read_chars = 0;

    for (int i = 0; i < size; ++i, offset += read_chars)
    {
        /* Пробелы между числами в строке останавливают чтение sscanf. Нужен отступ. */
        sscanf(input_numbers + offset, "%u%n", &number, &read_chars);
        if ((!i && !number) || !number)
        {
            break;
        }
        counter_even += number % DIVIDER ? 0 : 1;
    }

    return counter_even;
}