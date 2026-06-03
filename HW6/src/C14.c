/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C14: Составить функцию логическую функцию, которая определяет, 
 *      верно ли, что сумма его цифр – четное число.
 */
#include <stdio.h>
#include <stdlib.h>

const char* is_even_sum_digits_of_number (unsigned number)
{
    enum { TWO = 2, TEN = 10 };
    unsigned sum = 0;
    while (number)
    {
        sum += number % TEN;
        number /= TEN;
    }
    return sum % TWO ? "NO" : "YES";
}

#ifndef TEST_DEF_HW6
int main (void)
{
    unsigned input_num = 0;
    scanf("%u", &input_num);
    printf("%s\n", is_even_sum_digits_of_number(input_num));
    return EXIT_SUCCESS;
}
#endif