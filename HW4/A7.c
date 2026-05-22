/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A7: Ввести два числа и вывести их в порядке возрастания.
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
    int first = a, second = b;
    if (a > b)
    {
        first = b;
        second = a;
    } else if (b > a)
    {
        first = a;
        second = b;
    }
    printf("%d %d\n", first, second);

    return EXIT_SUCCESS;
}
#endif