/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B14: Дана последовательность ненулевых целых чисел, 
 *      в конце последовательности число 0. Посчитать количество чисел
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef unsigned uint;

#ifndef TEST_DEF_HW5
int main (void)
{
    int counter_numbers = -1;
    uint number = UINT32_MAX;
    while (number && scanf("%u", &number))
    {
        counter_numbers++;
    }
    printf("%u\n", counter_numbers);

    return EXIT_SUCCESS;
}
#endif