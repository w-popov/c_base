/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B16: Составить программу для вычисления НОД с помощью 
 *      алгоритма Евклида. Даны два натуральных числа. 
 *      Найти наибольший общий делитель.
 */
#include <stdio.h>
#include <stdlib.h>

/**
 * Н.О.Д алгоритм Евклида
 */
int euclid_algorithm (int number_a, int number_b)
{
    while (number_b)
    {
        int temporary = number_b;
        number_b = number_a % number_b;
        number_a = temporary;
    }
    return number_a > 0 ? number_a : -number_a;
}

#ifndef TEST_DEF_HW5
int main (void)
{
    int number_a = 0, number_b = 0;
    scanf("%d %d", &number_a, &number_b);
    printf("%d\n", euclid_algorithm(number_a, number_b));

    return EXIT_SUCCESS;
}
#endif