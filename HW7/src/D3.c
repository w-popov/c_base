/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D3: Дано целое не отрицательное число N.
 *      Выведите все его цифры по одной, в обратном порядке,
 *      разделяя их пробелами или новыми строками
 */
#include <stdio.h>
#include <stdlib.h>
// #include "HW7.h"

void print_digits_reverce_rec (int number)
{
    printf("%d ", number % 10);

    if (number / 10 != 0)
    {
        print_digits_reverce_rec(number / 10);
    }
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    print_digits_reverce_rec(input_number);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif