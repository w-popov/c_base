/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E11: Считать массив из 10 элементов и отсортировать 
 *      его по последней цифре (по возрастанию). 
 */
#include <stdio.h>
#include <stdlib.h>

int *fill_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        scanf("%d", &array[i]);
    }
    return array;
}

typedef int (*Compare)(int, int);

int compate_last_digit (int a, int b)
{
    return (a % 10) < (b % 10);
}

/**
 * Сортировка массива int[] вставками по последней цифре
 */
int *sort_X_array (int *array, const int size_array, Compare cmp)
{
    for (int i = 1; i < size_array; ++i)
    {
        int key = array[i];
        int j = i - 1;
        while (j >= 0 && cmp(key, array[j]))
        {
            array[j + 1] = array[j];
            j--;
        }
        array[j + 1] = key;
    }
    return array;
}

void print_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        printf("%d%s", array[i], i < size_array - 1 ? " " : "\n");
    }
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};

    print_array (
        sort_X_array (
            fill_array (array, SIZE_ARR),
            SIZE_ARR,
            compate_last_digit
        ), 
        SIZE_ARR
    );
    
    return EXIT_SUCCESS;
}
#endif