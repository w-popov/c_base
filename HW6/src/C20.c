/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C20: Проверить строку состоящую из скобок "(" и ")" на корректность.
 */
#include <stdio.h>
#include <stdlib.h>

const char *is_valid_brackets (const char *text, size_t size)
{
    const char open = '(';
    const char closed = ')';
    int valid_brackets = 0;
    for (size_t i = 0; i < size; ++i)
    {
        if (!i && text[i] == closed)
        {
            return "NO";
        }
        if (!text[i])
        {
            break;
        }
        if (text[i] == open)
        {
            ++valid_brackets;
        }
        else if (text[i] == closed)
        {
            --valid_brackets;
        }
        else
        {
            continue;
        }
    }
    return valid_brackets ? "NO" : "YES";
}

#ifndef TEST_DEF_HW6
int main (void)
{
    enum { SIZE_ARRAY_S = 1028 };
    char text[SIZE_ARRAY_S] = {'\0'};
    scanf(" %[^.]", text);
    printf("%s\n", is_valid_brackets((const char *)text, (size_t)SIZE_ARRAY_S));
    return EXIT_SUCCESS;
}
#endif