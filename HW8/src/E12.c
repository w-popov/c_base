/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E12: Считать массив из 10 элементов и отсортировать 
 *      первую половину по возрастанию, а вторую – по убыванию.
 */
#include <stdio.h>
#include <stdlib.h>

static int *fill_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        scanf("%d", &array[i]);
    }
    return array;
}

typedef int (*Compare)(int, int);

int compare_asc (int a, int b)
{
    return a < b;
}

int compare_desc (int a, int b)
{
    return b < a;
}

/**
 * Сортировка массива int[] вставками
 */
static int *sort_X_array (int *array, const int size_array, Compare cmp)
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

/**
 * Сортировка по условию задачи
 */
int *sort_array (int *array, const int size_array)
{
    const int half_size_array = size_array / 2;
    sort_X_array(array, half_size_array, compare_asc);
    sort_X_array(array + half_size_array, half_size_array, compare_desc);
    return array;
}

static void print_array (int *array, const int size_array)
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
        sort_array (
            fill_array (array, SIZE_ARR),
            SIZE_ARR
        ), 
        SIZE_ARR
    );
    
    return EXIT_SUCCESS;
}
#endif