/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B13: Посчитать количество четных и нечетных цифр числа
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef unsigned uint;

/**
 * Посчитать количество четных и нечетных цифр числа.
 * Результаты передаются по указателям *even, *odd через 
 * параметры в место вызова функции.
 */
int calculate_even_odd_digits_number (uint* even, uint* odd, uint _number)
{
    enum { ZERO, ONE, TWO, TEN=10 };

    if ( !(even && odd) )
        return EXIT_FAILURE;

    if (_number < TEN)
    {
        if ( _number % TWO )
        {
            *even = ZERO;
            *odd = ONE;
        } else {
            *even = ONE;
            *odd = ZERO;
        }
        return EXIT_SUCCESS;
    }

    for (uint number = _number; number; number /= TEN)
    {
        uint remainder = number % TEN;
        if ( remainder % TWO )
            ++(*odd);
        else
            ++(*even);
    }
    return EXIT_SUCCESS;
}


#ifndef TEST_DEF_HW5
int main (void)
{
    uint number = 0;
    uint even = 0, odd = 0;
    scanf("%u", &number);

    if ( calculate_even_odd_digits_number(&even, &odd, number) )
        return EXIT_FAILURE;

    printf("%u %u\n", even, odd);
    
    return EXIT_SUCCESS;
}
#endif