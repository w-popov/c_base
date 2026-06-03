/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C9: Составить функцию вычисления N!.
 *     Использовать ее при вычислении факториала int factorial(int n)
 */
#include <stdio.h>
#include <stdlib.h>

int factorial (int n)
{
    enum { ONE = 1 };
    if (!n)
    {
        return ONE;
    }
    int fact_N = ONE;

    for (int i = ONE; i <= n; ++i)
    {
        fact_N = fact_N * i;
    }

    return fact_N;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_num = 0;
    scanf("%d", &input_num);
    printf("%d\n", factorial(input_num));
    return EXIT_SUCCESS;
}
#endif