/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F14: Составить функцию которая возвращает сумму элементов 
 *      в заданном отрезке [from, to] для массива.
 *      Прототип функции int sum_between_ab(int from, int to, int size, int a[])
 */
#include <stdio.h>
#include <stdlib.h>

int sum_between_ab(int from, int to, int size, int a[])
{
    if (from > to)
    {
        int tmp = from;
        from = to;
        to = tmp;
    }
    int sum = 0;
    for (int i = 0; i < size; ++i)
    {
        sum += (a[i] >= from && a[i] <= to) ? a[i] : 0;
    }
    return sum;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif