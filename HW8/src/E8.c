/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E8: Считать массив из 12 элементов и выполнить инверсию для каждой трети массива
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

int *third_reverse_array (int *array, const int size_array)
{
    int third = (size_array / 3);
    int i_begin = 0, i_end = third - 1;
    int k_begin = third, k_end = (third * 2) - 1;
    int m_begin = third * 2, m_end = size_array - 1;
    for (; i_begin < third / 2;)
    {
        int bucket_i = array[i_begin];
        array[i_begin++] = array[i_end];
        array[i_end--] = bucket_i;

        int bucket_k = array[k_begin];
        array[k_begin++] = array[k_end];
        array[k_end--] = bucket_k;

        int bucket_m = array[m_begin];
        array[m_begin++] = array[m_end];
        array[m_end--] = bucket_m;
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
    enum { SIZE_ARR = 12 };
    int array[SIZE_ARR] = {0};

    print_array (
        third_reverse_array (
            fill_array (array, SIZE_ARR),
            SIZE_ARR), 
        SIZE_ARR);
    
    return EXIT_SUCCESS;
}
#endif