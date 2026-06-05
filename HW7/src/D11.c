/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D11: Дано натуральное число N. Посчитать количество «1» в двоичной записи
 * числа.
 */
#include <stdio.h>
#include <stdlib.h>

int how_many_truths (int num)
{
    if (!num)
    {
        return 0;
    }
    return (num % 2) + how_many_truths(num / 2);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%d\n", how_many_truths(input_number));
    return EXIT_SUCCESS;
}
#endif