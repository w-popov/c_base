/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E3: Считать массив из 10 элементов и найти в нем
 *     максимальный и минимальный элементы и их номера.
 * out: 4 целых числа через пробел:
 *      номер максимума, максимум, номер минимума, минимум
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

struct MinMax
{
    unsigned max_index;
    int max;
    unsigned min_index;
    int min;
};

static int* fill_array (int *array, const int size_array)
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

struct MinMax min_max_from_array (int *array, const int size)
{
    int min = INT32_MAX, max = INT32_MIN;
    unsigned min_i = 0, max_i = 0;
    for (int i = 0; i < size; ++i)
    {
        if (array[i] < min)
        {
            min = array[i];
            min_i = i + 1;
        }
        if (array[i] > max)
        {
            max = array[i];
            max_i = i + 1;
        }
    }
    return (struct MinMax){max_i, max, min_i, min};
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    struct MinMax minmax = min_max_from_array(
                            fill_array(array, SIZE_ARR), 
                            SIZE_ARR);
        
    printf("%u %d %u %d\n", 
        minmax.max_index, 
        minmax.max, 
        minmax.min_index,
        minmax.min
    );

    return EXIT_SUCCESS;
}
#endif