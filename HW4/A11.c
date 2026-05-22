/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A11: Напечатать сумму максимума и минимума.
 */

#include <stdio.h>
#include <stdlib.h>

/**
 * Вернуть сумму максимума и минимума из 5 чисел.
 * (чтобы определить здесь функции из предыдущ. заданий,
 * то нужно их выносить в отдельный файл. линковщик тестов
 * ругается)
 */
int sum_min_max_from_5_numbers(int a, int b, int c, int d, int e)
{
    int max = a, min = a;
    if (b > max) max = b;
    if (c > max) max = c;
    if (d > max) max = d;
    if (e > max) max = e;
    
    if (b < min) min = b;
    if (c < min) min = c;
    if (d < min) min = d;
    if (e < min) min = e;
    
    return min + max;
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
    printf("%d", sum_min_max_from_5_numbers(a, b, c, d, e));

    return EXIT_SUCCESS;
}
#endif