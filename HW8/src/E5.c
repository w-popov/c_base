/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E5: Считать массив из 10 элементов и посчитать 
 *     сумму положительных элементов массива 
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

int sum_positive_array (int *array, const int size_array)
{
    int sum = 0;
    for (int i = 0; i < size_array; ++i)
    {
        sum += array[i] > 0 ? array[i] : 0;
    }
    return sum;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    printf("%d\n",
           sum_positive_array(
            fill_array(array, SIZE_ARR), 
            SIZE_ARR)
        );
    return EXIT_SUCCESS;
}
#endif