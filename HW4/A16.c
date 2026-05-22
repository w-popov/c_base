/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A16: Ввести три числа и определить, верно ли, что 
 * они вводились в порядке возрастания.
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
    const char *yes = "YES", *no = "NO";
    int a = 0, b = 0, c = 0;
    scanf("%d %d %d", &a, &b, &c); 
    printf("%s\n", (a < b) && (b < c) ? yes : no);

    return EXIT_SUCCESS;
}
#endif