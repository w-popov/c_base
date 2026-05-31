/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B1: Ввести натуральное число вывести квадраты и кубы 
 *     всех чисел от 1 до этого числа. Число не превосходит 100.
 */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#ifndef TEST_DEF_HW5
int main (void)
{
    uint16_t enter_number = 0;
    scanf("%"SCNu16, &enter_number);

    for (uint16_t i = 1; i <= enter_number; ++i)
    {
        printf("%" PRIu32 " " "%" PRIu32 " " "%" PRIu32 "\n", i, (i * i), (i * i * i));
    }

    return EXIT_SUCCESS;
}
#endif