/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D1: Составить рекурсивную функцию, печать всех чисел от 1 до N
 */
#include <stdio.h>
#include <stdlib.h>
//  #include "HW7.h"

void print_range_rec (int number)
{
    if (!number)
    {
        return;
    }
    print_range_rec(number - 1);
    printf("%d ", number);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    print_range_rec(input_number);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif