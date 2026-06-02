/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C8: Составить функцию, которая переводит латинскую
 *     строчную букву в заглавную
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/**
 * Перевести латинскую строчную букву в заглавную
 */
const char *to_uppercase (char *text, char *upper, const int str_size)
{
    memset(upper, '\0', str_size);
    enum { UPPER_SHIFT = 5 };

    uint8_t mask = 1 << UPPER_SHIFT;
    for (int i = 0; i < str_size; ++i)
    {
        if (!text[i])
        {
            break;
        }
        if (text[i] >= 97 && text[i] <= 122)
        {
            upper[i] = text[i] & ~mask;
            continue;
        }
        upper[i] = text[i];
    }
    return (const char *)upper;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    enum { STR_SIZE = 512 };
    char string[STR_SIZE] = {'\0'};
    char upper[STR_SIZE];
    scanf(" %[^.]", string);
    printf("%s\n", to_uppercase(string, upper, STR_SIZE));
    return EXIT_SUCCESS;
}
#endif