/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D6: Дана строка заканчивающаяся символом точка '.'
 *      Напечатать строку наоборот. Реализуйте рекурсивную функцию,
 *      которая считывает и печатает строку наоборот. void reverse_string()
 */
#include <stdio.h>
#include <stdlib.h>

void reverse_string ()
{
    char c = '\0';
    scanf("%c", &c);
    if (c == '.')
    {
        return;
    }
    reverse_string();
    printf("%c", c);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    reverse_string();
    printf("\n");
    return EXIT_SUCCESS;
}
#endif