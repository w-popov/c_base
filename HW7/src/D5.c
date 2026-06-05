/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D5: Перевести натуральное число в двоичную систему счисления.
 *     Необходимо реализовать рекурсивную функцию.
 */
#include <stdio.h>
#include <stdlib.h>

unsigned dec_to_bin_rec (unsigned num)
{
    if (!num)
    {
        return 0;
    }
    unsigned rest = num % 2;
    return dec_to_bin_rec(num / 2) * 10 + rest;
}

#ifndef TEST_DEF_HW7
int main (void)
{
    unsigned input_number = 0;
    scanf("%u", &input_number);
    printf("%u\n", dec_to_bin_rec(input_number));
    return EXIT_SUCCESS;
}
#endif