/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B4: Ввести целое signed число и определить,
 *     верно ли, что в нём ровно 3 цифры
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW5.h"

#ifndef TEST_DEF_HW5
int main (void)
{
    unsigned number = 0;
    scanf("%u", &number);
    printf("%s\n", is_it_number_have_3_digit(number));

    return EXIT_SUCCESS;
}
#endif

const char *is_it_number_have_3_digit (unsigned number)
{
    unsigned counter_digit = 0;
    enum { III = 3 };

    while (number)
    {
        number /= 10;
        ++counter_digit;
    }

    return counter_digit == III ? "YES" : "NO";
}