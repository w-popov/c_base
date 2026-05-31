/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B8: Ввести целое число и определить, верно ли, 
 *     что в нём ровно одна цифра «9»
 */

#include <stdio.h>
#include <stdlib.h>

const char* is_nine_digit_in_number (int number)
{
    enum { SIZE = 16 };
    char str_array_number[SIZE] = {'\0'};
    snprintf(str_array_number, SIZE, "%d", number);
    int nine_counter = 0;
    for (int i = 0; i < SIZE; ++i)
    {
        nine_counter += ( (str_array_number[i] - '0') == 9 ) ? 1 : 0;
    }
    return (nine_counter != 1) ? "NO" : "YES";
}

#ifndef TEST_DEF_HW5
int main (void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%s\n", is_nine_digit_in_number(number));

    return EXIT_SUCCESS;
}
#endif