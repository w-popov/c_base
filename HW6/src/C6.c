/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C6: Когда создатель шахмат, древнеиндийский мудрец и
 *     математик Сисса бен Дахир ...
 */
#include <stdio.h>
#include <stdlib.h>

unsigned long long how_many_grains (unsigned cell)
{
    if (cell > 0 && cell < 3)
    {
        return (unsigned long long)cell;
    }
    enum { ONE = 1, TWO = 2 };
    unsigned long long grains = ONE;
    for (unsigned i = grains; i < cell; ++i)
    {
        grains *= TWO;
    }
    return grains;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    unsigned cell = 0;
    scanf("%u", &cell);
    printf("%llu\n", how_many_grains(cell));
    return EXIT_SUCCESS;
}
#endif