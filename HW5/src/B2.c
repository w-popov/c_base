/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B2: Ввести два целых числа a и b (a ≤ b) и вывести квадраты всех чисел от a до b.
 */

#include <stdio.h>
#include <stdlib.h>

int square_power (int number) { return number * number; }

#ifndef TEST_DEF_HW5
int main (void)
{
    int numberA = 0, numberB = 0;
    scanf("%d %d", &numberA, &numberB);

    for (int i = numberA; i <= numberB; ++i)
    {
        printf("%d%s" , square_power(i), i < numberB ? " " : "\n");
    }

    return EXIT_SUCCESS;
}
#endif