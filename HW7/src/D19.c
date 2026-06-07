/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D19: Дана строка из английских символов, пробелов и
 *      знаков препинания. В конце строки символ точка.
 *      Необходимо реализовать рекурсивную функцию, которая
 *      считывает данную строку со стандартного потока ввода и
 *      возвращает целое число - количество символов 'a' в данной строке.
 *      int acounter(void)
 */
#include <stdio.h>
#include <stdlib.h>

int acounter (void)
{
    char ch = '\0';
    int count_a = 0;
    scanf("%c", &ch);
    if (ch == '.')
    {
        return 0;
    }
    if (ch == 'a')
    {
        ++count_a;
    }
    return count_a + acounter();
}

#ifndef TEST_DEF_HW7
int main (void)
{
    printf("%d\n", acounter());
    return EXIT_SUCCESS;
}
#endif