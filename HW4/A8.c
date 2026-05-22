/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A8: Ввести три числа и найти наибольшее из них
 */

#include <stdio.h>
#include <stdlib.h>

/**
 * Вернуть наибольшее из трёх чисел a, b и c.
 */
int find_max(int a, int b, int c)
{
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

/**
 * TEST_DEF_HW определяется в ./tests/Makefile для компиляции тестов. 
 * Если он определён, то тесты смогут использовать функции из этого файла.
 */
#ifndef TEST_DEF_HW
int main(void)
{
    int a = 0, b = 0, c = 0;
    scanf("%d %d %d", &a, &b, &c);
    printf("%d\n", find_max(a, b, c));
    
    return EXIT_SUCCESS;
}
#endif