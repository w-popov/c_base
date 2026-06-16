/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E7: Считать массив из 10 элементов и выполнить инверсию 
 *     отдельно для 1-ой и 2-ой половин массива. Необходимо 
 *     изменить считанный массив и напечатать его одним циклом. 
 */
#include <stdio.h>
#include <stdlib.h>

static int *fill_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        if (scanf("%d", &array[i]) != 1)
        {
            printf("Error scanf\n");
            exit(EXIT_FAILURE);
        }
    }
    return array;
}

int *half_reverse_array (int *array, const int size_array)
{
    int half = (size_array / 2);
    int i_begin = 0, i_end = half - 1;
    int k_begin = half, k_end = size_array - 1;
    for (; i_begin < (half / 2); ++i_begin, --i_end, ++k_begin, --k_end)
    {
        int bucket_i = array[i_begin];
        array[i_begin] = array[i_end];
        array[i_end] = bucket_i;
        int bucket_k = array[k_begin];
        array[k_begin] = array[k_end];
        array[k_end] = bucket_k;
    }
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
        half_reverse_array (
            fill_array (array, SIZE_ARR),
            SIZE_ARR), 
        SIZE_ARR);
    
    return EXIT_SUCCESS;
}
#endif