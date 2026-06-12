/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E9: Считать массив из 10 элементов и выполнить циклический сдвиг ВПРАВО на 1 
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

// 1 2 3 4 5 6 7 8 9 10
int *rshift_array (int *array, int shift, const int size)
{
    shift = ((shift % size) + size) % size;
    int shift_array[size];
    for (int i = 0; i < size; i++) 
    {
        shift_array[(i + shift) % size] = array[i];
    }
    for (int k = 0; k < size; ++k)
    {
        array[k] = shift_array[k];
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
        rshift_array (
            fill_array (&array[0], SIZE_ARR),
            1,
            SIZE_ARR), 
        SIZE_ARR);
    
    return EXIT_SUCCESS;
}
#endif