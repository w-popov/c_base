/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F13: Составить функцию которая возвращает количество 
 *      элементов на заданном отрезке [from, to] для массива.
 *      Прототип функции int count_between(int from, int to, int size, int a[]) 
 */
#include <stdio.h>
#include <stdlib.h>

int count_between(int from, int to, int size, int a[])
{
    int count = 0;
    for (int i = 0; i < size; ++i)
    {
        count += (a[i] >= from && a[i] <= to) ? 1 : 0;
    }
    return count;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif