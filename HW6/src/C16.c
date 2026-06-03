/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C16: Составить функцию логическую функцию, которая
 *      определяет, верно ли, что число простое.
 *      int is_prime(int n)
 */
#include <stdio.h>
#include <stdlib.h>

int is_prime (int n)
{
    enum { ZERO, ONE, TWO };
    if (n < TWO)
    {
        return ZERO;
    }
    for (int counter = TWO; (counter * counter) <= n; ++counter)
    {
        if (!(n % counter))
        {
            return ZERO;
        }
    }
    return ONE;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_num = 0;
    scanf("%d", &input_num);
    printf("%s\n", is_prime(input_num) ? "YES" : "NO");
    return EXIT_SUCCESS;
}
#endif