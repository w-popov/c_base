/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A10: Ввести пять чисел и найти наименьшее из них.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/**
 * Вернуть наименьшее из 5 чисел
 */
int min_from_5_numbers(int a, int b, int c, int d, int e)
{
    int min = a;
    if (b < min) min = b;
    if (c < min) min = c;
    if (d < min) min = d;
    if (e < min) min = e;
    
    return min;
}

/**
 * TEST_DEF_HW определяется в Makefile для компиляции тестов в ./tests. 
 * Если он определён, то main в исходниках си не будет компилироваться, 
 * и тесты смогут использовать функции из этих исходников.
 */
#ifndef TEST_DEF_HW
int main(void)
{
    int a = 0, b = 0, c = 0, d = 0, e = 0;
    scanf("%d %d %d %d %d", &a, &b, &c, &d, &e);
    printf("%d", min_from_5_numbers(a, b, c, d, e));

    return EXIT_SUCCESS;
}
#endif