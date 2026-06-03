/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C19: Составить функцию, которая преобразует текущий
 *      символ цифры в число. Написать программу считающую
 *      сумму цифр в тексте. int digit_to_num(char c)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int digit_to_num (char c)
{
    enum { START_ASCII_NUM = 48, END_ASCII_NUM = 57 };
    return c >= START_ASCII_NUM && c <= END_ASCII_NUM ? c - '0' : NAN;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    enum { SIZE_ARRAY_S = 512 };
    char text[SIZE_ARRAY_S] = {'\0'};
    int sum_nums = 0;
    scanf(" %[^.]", text);
    for (char *pt_c = text; *pt_c; ++pt_c)
    {
        sum_nums += isnan((float)digit_to_num(*pt_c)) ? 0 : digit_to_num(*pt_c);
    }
    printf("%d\n", sum_nums);
    return EXIT_SUCCESS;
}
#endif