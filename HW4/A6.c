/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A6: Ввести два числа и найти их разность.
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
    int a = 0, b = 0;
    scanf("%d %d", &a, &b);
    printf("%d\n", a - b);

    return EXIT_SUCCESS;
}
#endif