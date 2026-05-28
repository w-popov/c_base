/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B11: Ввести целое число и «перевернуть» его, так 
 *      чтобы первая цифра стала последней и т.д
 */

#include <stdio.h>
#include <stdlib.h>

unsigned reverse_number(unsigned number)
{
    enum { DEC = 10 };
    if (number < DEC)
        return number;
        
    unsigned reversed_number = 0;
    while (number)
    {
        reversed_number = (reversed_number * DEC) + (number % DEC);
        number /= DEC;
    }

    return reversed_number;
}


#ifndef TEST_DEF_HW5
int main(void)
{
    unsigned number = 0;
    scanf("%u", &number);
    printf("%u\n", reverse_number(number));

    return EXIT_SUCCESS;
}
#endif