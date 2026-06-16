/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E13: Считать массив из 10 элементов и отобрать в другой
 *      массив все числа, у которых вторая с конца цифра
 *      (число десятков) – ноль.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

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

int *round_number_array (int *array, int *result_aray, const int size_array)
{
    for (int i = 0, j = 0; i < size_array; ++i)
    {
        result_aray[i] = INT32_MIN;
        if ( !((array[i] / 10) % 10) )
        {
            result_aray[j++] = array[i];
        }
    }
    return result_aray;
}

static void print_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        if (array[i] != INT32_MIN)
        {
            printf("%d ", array[i]);
        }
    }
    printf("\n");
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    int result_array[SIZE_ARR] = {0};

    print_array(
        round_number_array(
            fill_array(array, SIZE_ARR), 
            result_array, SIZE_ARR
        ),
        SIZE_ARR);

    return EXIT_SUCCESS;
}
#endif