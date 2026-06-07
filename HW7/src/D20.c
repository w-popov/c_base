/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D20: Написать рекурсивную функцию возведения целого числа n в степень p.
 *      int recurs_power(int n, int p)
 */
#include <stdio.h>
#include <stdlib.h>

int recurs_power (int n, int p)
{
    if (!p)
    {
        return 1;
    }
    return n * recurs_power(n, p - 1);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number_n = 0, input_number_p = 0;
    scanf("%d %d", &input_number_n, &input_number_p);
    printf("%d\n", recurs_power(input_number_n, input_number_p));
    return EXIT_SUCCESS;
}
#endif