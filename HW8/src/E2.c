/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E2: Ввести c клавиатуры массив из 5 элементов, найти минимальный из них. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int* fill_array (int *array, const int size_array)
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

int min_num_array (int *array, const int size_array)
{
    int min = INT32_MAX;
    for (int i = 0; i < size_array; ++i)
    {
        min = min < array[i] ? min : array[i];
    }
    return min;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 5 };
    int array[SIZE_ARR] = {0};
    printf("%d\n", min_num_array(fill_array(array, SIZE_ARR), SIZE_ARR));
    return EXIT_SUCCESS;
}
#endif