/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F4: Написать только одну функцию. Всю программу отправлять не надо. 
 *     Вывести в порядке возрастания цифры, входящие в строку. 
 *     Цифра - количество. Функция должно строго соответствовать прототипу:
 *     void print_digit(char s[])
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void print_digit(char s[])
{
    enum { MAX10 = 10 };
    size_t digits[MAX10] = {0};
    for (size_t i = 0; i < strlen(s); ++i)
    {
        if (isdigit(s[i]))
        {
            digits[s[i] - '0']++;
        }
    }
    for (size_t i = 0; i < MAX10; ++i)
    {
        if (digits[i])
        {
            printf("%zu %zu", i, digits[i]);
            printf("%s", i < MAX10 ? "\n" : "");
        }
    }
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { MAXLEN = 1028 };
    char str[MAXLEN] = {'\0'};
    if (scanf(" %[^\n]", &str[0]) != 1)
    {
        return EXIT_FAILURE;
    }
    print_digit(str);
    return EXIT_SUCCESS;
}
#endif