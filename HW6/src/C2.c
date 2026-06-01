/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C2: Составить функцию, возведение числа N в степень P.
 *     int power(n, p) и привести пример ее использования
 */
#include <stdio.h>
#include <stdlib.h>

int power (int n, unsigned p);

#ifndef TEST_DEF_HW6
int main (void)
{
    int number = 0;
    unsigned pow = 0;
    scanf("%d"
          "%u",
          &number, &pow);
    printf("%d\n", power(number, pow));
    return EXIT_SUCCESS;
}
#endif

int power (int n, unsigned p)
{
    int powered = 1;
    for (unsigned i = p; i; --i)
    {
        powered = powered * n;
    }
    return powered;
}