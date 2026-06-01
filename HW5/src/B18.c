/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B18: Вывести на экран ряд чисел Фибоначчи, состоящий из n элементов
 */
#include <stdio.h>
#include <stdlib.h>

#ifndef TEST_DEF_HW5
int main (void)
{
    unsigned number = 0;
    unsigned prev_a = 0, prev_b = 1, it_number = 1;
    scanf("%u", &number);

    if (number == 1)
    {
        printf("%u\n", number);
        return EXIT_SUCCESS;
    }

    for (unsigned i = 1; i <= number; ++i)
    {
        printf("%u ", it_number);
        it_number = prev_a + prev_b;
        prev_a = prev_b;
        prev_b = it_number;
    }
    printf("\n");

    return EXIT_SUCCESS;
}
#endif
