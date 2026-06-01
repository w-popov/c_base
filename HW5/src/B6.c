/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B6: Ввести целое число и определить, верно ли,
 * что в его записи есть две одинаковые цифры, стоящие рядом
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW5.h"

#ifndef TEST_DEF_HW5
int main (void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%s\n", is_both_equal_digits_in_number(number));

    return EXIT_SUCCESS;
}
#endif

const char *is_both_equal_digits_in_number (int number)
{
    int prev_digits = 0;
    if (number < 10)
    {
        return "NO";
    }

    do
    {
        prev_digits = (number % 10);
        number /= 10;
        if ((number % 10) == prev_digits)
        {
            return "YES";
        }

    } while (number);

    return "NO";
}