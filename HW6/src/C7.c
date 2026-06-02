/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C7: Составить функцию, которая переводит число N
 *     из десятичной системы счисления в P-ичную систему счисления.
 */
#include <stdio.h>
#include <stdlib.h>

/**
 * Перевод <number> из десятичной системы 
 * счисления в <p>-ичную систему счисления
 */
int p_notation (int number, unsigned p)
{
    enum { DEC = 10 };
    int result_number = 0, factor = 1;
    while (number)
    {
        int rest = number % p;
        result_number += (rest * factor);
        factor *= DEC;
        number /= p;
    }
    return result_number;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_num = 0;
    unsigned p = 0;
    scanf("%d %u", &input_num, &p);
    printf("%d\n", p_notation(input_num, p));
    return EXIT_SUCCESS;
}
#endif