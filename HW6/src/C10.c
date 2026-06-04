/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C10: Составить функцию, печать всех простых множителей числа.
 *      Использовать ее для печати всех простых множителей числа.
 *      void print_simple(int n)
 */
#include <stdio.h>
#include <stdlib.h>

void print_simple (int n)
{
    enum { TWO = 2, THREE = 3 };
    while ( !(n % TWO) )
    {
        printf("%d ", TWO);
        n /= TWO;
    }

    for (int i = THREE; (i * i) <= n; i += TWO)
    {
        while ( !(n % i) )
        {
            printf("%d ", i);
            n /= i;
        }
    }
    n > TWO ? printf("%d\n", n) : printf("\n");
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_num = 0;
    scanf("%d", &input_num);
    print_simple(input_num);
    return EXIT_SUCCESS;
}
#endif