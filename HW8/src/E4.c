/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E4: Считать массив из 10 элементов и найти в нем
 *     два максимальных элемента и напечатать их сумму
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int *fill_array (int *array, const int size_array)
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

int sum_two_max_from_array (int *array, const int size_array)
{
    int max_1 = INT32_MIN, max_2 = INT32_MIN;
    if (array[0] >= array[1])
    {
        max_1 = array[0];
        max_2 = array[1];
    }
    else
    {
        max_1 = array[1];
        max_2 = array[0];
    }
    for (int i = 2; i < size_array; ++i)
    {
        if (array[i] > max_1)
        {
            max_2 = max_1;
            max_1 = array[i];
        }
        else if (array[i] > max_2)
        {
            max_2 = array[i];
        }
    }
    return max_1 + max_2;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    printf("%d\n",
           sum_two_max_from_array(
            fill_array(array, SIZE_ARR), 
            SIZE_ARR)
        );
    return EXIT_SUCCESS;
}
#endif