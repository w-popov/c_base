/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C18: Составить логическую функцию, которая определяет,
 *      что текущий символ это цифра. Написать программу
 *      считающую количество цифр в тексте. int is_digit(char c)
 */
#include <stdio.h>
#include <stdlib.h>

int is_digit (char c)
{
    enum { NO, YES, START_ASCII_NUM = 48, END_ASCII_NUM = 57 };
    return c >= START_ASCII_NUM && c <= END_ASCII_NUM ? YES : NO;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    enum { SIZE_ARRAY_S = 512 };
    char text[SIZE_ARRAY_S] = {'\0'};
    int count_nums = 0;
    scanf(" %[^.]", text);
    for (char* pt_c = text; *pt_c; ++pt_c)
    {
        count_nums += is_digit(*pt_c) ? 1 : 0;
    }
    printf("%d\n", count_nums);
    return EXIT_SUCCESS;
}
#endif