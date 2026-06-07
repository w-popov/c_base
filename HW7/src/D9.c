/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D9: Дано натуральное число N. Вычислите сумму его цифр.
 *     Необходимо составить рекурсивную функцию. int sum_digits(int n)
 */
#include <stdio.h>
#include <stdlib.h>

int sum_digits (int n)
{
    if (!n)
    {
        return 0;
    }
    return n % 10 + sum_digits(n / 10);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%d\n", sum_digits(input_number));
    return EXIT_SUCCESS;
}
#endif