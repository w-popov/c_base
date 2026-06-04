/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D2: Составить рекурсивную функцию, которая определяет
 *     сумму всех чисел от 1 до N
 */
#include <stdio.h>
#include <stdlib.h>
// #include "HW7.h"

int sum_range_numbers_rec (int number)
{
    if (!number)
    {
        return 0;
    }
    
    return number + sum_range_numbers_rec(number - 1);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%u\n", sum_range_numbers_rec(input_number));
    return EXIT_SUCCESS;
}
#endif