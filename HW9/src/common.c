#include "HW9.h"
#include <stdio.h>
#include <stdlib.h>
 
/**
 * F1: Функция сортирует массив по возрастанию.
 */
void sort_array(int size, int a[])
{
    for (int i = 1; i < size; ++i)
    {
        int key = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > key)
        {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = key;
    }
}

/**
 * F2: Функция ставит в начало массива все четные элементы, а в конец – все нечетные
 */
void sort_even_odd(int n, int a[])
{
    int *array = calloc(n, sizeof(int));
    int *array_even = array;
    
    for (int i = 0; i < n; ++i)
    {
        if (a[i] % 2 == 0)
        {
            *array_even++ = a[i];
        }
    }
    for (int i = 0; i < n; ++i)
    {
        if (a[i] % 2 != 0)
        {
            *array_even++ = a[i];
        }
    }
    for (int j = 0; j < n; ++j)
    {
        a[j] = array[j];
    }
    free(array);
}