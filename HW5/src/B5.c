/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B5: Ввести целое число и найти сумму его цифр.
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW5.h"

#ifndef TEST_DEF_HW5
int main (void)
{
    unsigned number = 0;
    scanf("%u", &number);
    printf("%u\n", summ_digits_of_number(number));
    
    return EXIT_SUCCESS;
}
#endif


unsigned summ_digits_of_number (unsigned number)
{
    unsigned summ_digits = 0;
    while (number)
    {
        summ_digits += (number % 10);
        number /= 10;
    }
    return summ_digits;
}