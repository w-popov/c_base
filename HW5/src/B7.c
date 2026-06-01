/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B7: Ввести целое число и определить, верно ли, что в его
 *     записи есть  две одинаковые цифры, НЕ обязательно стоящие рядом
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW5.h"

#ifndef TEST_DEF_HW5
int main (void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%s\n", is_equal_digits_in_number(number));

    return EXIT_SUCCESS;
}
#endif

const char *is_equal_digits_in_number (int number)
{
    if (number < 10)
    {
        return "NO";
    }
    enum { SIZE = 16 };
    char str_array_number[SIZE] = {'\0'};
    int int_array_number[SIZE] = {0};
    snprintf(str_array_number, SIZE, "%d", number);

    for (int i = 0; i < SIZE; ++i)
    {
        int_array_number[(str_array_number[i] - '0')] += 1;
    }
    for (int i = 0; i < SIZE; ++i)
    {
        if (int_array_number[i] > 1)
        {
            return "YES";
        }
    }

    return "NO";
}