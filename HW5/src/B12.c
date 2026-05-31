/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B12: Ввести натурального числа с клавиатуры. 
 *      Программа должна определить наименьшую и наибольшую цифры, 
 *      которые входят в состав данного натурального числа
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "HW5.h"

typedef unsigned uint;

#ifndef TEST_DEF_HW5
int main (void)
{
    uint number = 0;
    uint min_digit = UINT32_MAX, max_digit = 0;
    scanf("%u", &number);
    
    if ( calculate_min_max_digits_number(&min_digit, &max_digit, number) )
    return EXIT_FAILURE;
    
    printf("%u %u\n", min_digit, max_digit);
    
    return EXIT_SUCCESS;
}
#endif


/**
 * Найти минимальную и максимальную цифры числа _number.
 * Результаты передаются по указателям *min, *max через 
 * параметры в место вызова функции.
 */
int calculate_min_max_digits_number (uint* min, uint* max, uint _number)
{
    enum { TEN=10 };

    if ( !(min && max) )
        return EXIT_FAILURE;

    if (_number < 10)
    {
        *min = _number;
        *max = _number;
        return EXIT_SUCCESS;
    }
    for (uint number = _number; number; number /= TEN)
    {
        uint remainder = number % TEN;
        if ( remainder < *min )
            *min = remainder;
        if ( remainder > *max)
            *max = remainder;
    }
    return EXIT_SUCCESS;
}