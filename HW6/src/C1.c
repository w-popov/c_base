/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C1: Составить функцию, модуль числа и привести пример ее использования.
 */
#include <stdio.h>
#include <stdlib.h>

unsigned mod (int number);

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%u\n", mod(input_number));
    return EXIT_SUCCESS;
}
#endif

unsigned mod (int number)
{
    return number < 0 ? (unsigned)(-number) : (unsigned)number;
}