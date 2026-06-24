#include "HW9.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
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

/**
 * F8: Определить пропущенное число в последовательности
 */
int missing_number (int *array, const int size)
{
    int *cp_array = malloc(size * sizeof(int));
    if (cp_array == NULL)
    {
        perror("Error alloc mem in missing_number(int[], const int)\n");
        return -1;
    }
    /* Копия массива */
    memcpy(cp_array, array, size * sizeof(int));

    /* Сортировка вставками */
    for (int i = 1; i < size; ++i)
    {
        int key = cp_array[i];
        int j = i - 1;
        while (j >= 0 && cp_array[j] > key)
        {
            cp_array[j + 1] = cp_array[j];
            j--;
        }
        cp_array[j + 1] = key;
    }
    /* Поиск */
    for (int i = 1; i < size; ++i)
    {
        if (cp_array[i] != (cp_array[i-1] + 1))
        {
            return cp_array[i] - 1;
        }
    }
    
    free(cp_array);
    return 0;
}