/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E14: Считать массив из 10 элементов и выделить в 
 *     другой массив все числа, которые встречаются более
 *     одного раза. В результирующем массиве эти числа 
 *     не должны повторяться. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int *fill_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        scanf("%d", &array[i]);
    }
    return array;
}

int *unique_array (int *array, int *uniq_arr, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        uniq_arr[array[i]] += 1;
    }
    return uniq_arr;
}

void print_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        if (array[i] > 1)
        {
            printf("%d ", i);
        }
    }
    printf("\n");
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10, UNIQ_ARR_SIZE = INT16_MAX };
    int array[SIZE_ARR] = {0};
    int result_array[UNIQ_ARR_SIZE] = {0};
    
    print_array(
        unique_array(
            fill_array(array, SIZE_ARR), 
            result_array, 
            SIZE_ARR
        ),
        UNIQ_ARR_SIZE
    );
    
    return EXIT_SUCCESS;
}
#endif