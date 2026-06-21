/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F15: Составить функцию которая определяет в массиве, 
 *      состоящем из положительных и отрицательных чисел, 
 *      сколько элементов превосходят по модулю максимальный элемент.
 *      Прототип функции int count_bigger_abs(int n, int a[]);
 */
#include <stdio.h>
#include <stdlib.h>

int count_bigger_abs(int n, int a[])
{
    int max_element = -2147483648;
    int counter = 0;
    for (int i = 0; i < n; ++i)
    {
        max_element = a[i] > max_element ? a[i] : max_element;
    }
    for (int i = 0; i < n; ++i)
    {
        int value = a[i] < 0 ? -a[i] : a[i];
        counter += value > max_element ? 1 : 0;
    }
    return counter;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif