/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F12: Составить функцию которая меняет в массиве минимальный и 
 *      максимальный элемент местами.
 *      Прототип функции void change_max_min(int size, int a[])
 */
#include <stdio.h>
#include <stdlib.h>

void change_max_min(int size, int a[])
{
    int min = 0, max = 1, temp_swap = 0;;
    for (int i = 0; i < size; ++i)
    {
        min = a[i] < a[min] ? i : min;
        max = a[i] > a[max] ? i : max; 
    }
    temp_swap = a[min];
    a[min] = a[max];
    a[max] = temp_swap;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif