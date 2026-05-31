/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B3: Ввести два целых числа a и b (a ≤ b) и вывести
 *     сумму квадратов всех чисел от a до b.
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW5.h"

#ifndef TEST_DEF_HW5
int main (void)
{
    int a = 0, b = 0;
    scanf("%d %d", &a, &b);
    printf("%d\n", summ_squares_range(a, b));
    
    return EXIT_SUCCESS;
}
#endif


int summ_squares_range (int a, int b)
{
    int summ_squares = 0;
    for (int i = a; i <= b; ++i)
    {
        summ_squares += (i * i);
    }
    return summ_squares;
}
