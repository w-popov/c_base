/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B9: Ввести целое число и определить, верно ли, что все его цифры четные
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW5.h"

#ifndef TEST_DEF_HW5
int main (void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%s\n", is_even_all_digits_in_number(number));
    
    return EXIT_SUCCESS;
}
#endif


const char* is_even_all_digits_in_number (int number)
{
    enum { SIZE = 16 };
    char str_array_number[SIZE] = {'\0'};
    snprintf(str_array_number, SIZE, "%d", number);
    for (int i = 0; i < SIZE; ++i)
    {
        if ( (str_array_number[i] - '0') % 2 != 0 )
            return "NO";
    }
    return "YES";
}