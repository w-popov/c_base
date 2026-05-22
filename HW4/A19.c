/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A19: Даны стороны треугольника a, b, c. 
 * Определить существует ли такой треугольник.
 */

#include <stdio.h>
#include <stdlib.h>

/**
 * TEST_DEF_HW определяется в Makefile для компиляции тестов в ./tests. 
 * Если он определён, то main в исходниках си не будет компилироваться, 
 * и тесты смогут использовать функции из этих исходников.
 */
#ifndef TEST_DEF_HW
int main(void)
{
    int a = 0, b = 0, c = 0;
    _Bool is_triangle = 0;
    const char *yes = "YES", *no = "NO";
    scanf("%d %d %d", &a, &b, &c);

    is_triangle =   ((a > 0) && (b > 0) && (c > 0)) &&
                    ((a + b) > c) && 
                    ((a + c) > b) && 
                    ((b + c) > a);
    printf("%s\n", is_triangle ? yes : no);

    return EXIT_SUCCESS;
}
#endif