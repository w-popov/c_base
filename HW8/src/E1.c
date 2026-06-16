/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E1: Ввести c клавиатуры массив из 5 элементов,
 *     найти среднее арифметическое всех элементов массива.
 */
#include <stdio.h>
#include <stdlib.h>

static float* fill_array (float *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        if (scanf("%f", &array[i]) != 1)
        {
            printf("Error scanf\n");
            exit(EXIT_FAILURE);
        }
    }
    return array;
}

float average_array (float *array, const int size_array)
{
    float average = 0.0f;
    for (int i = 0; i < size_array; ++i)
    {
        average += array[i];
    }
    return average / size_array;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 5 };
    float array[SIZE_ARR] = {0.0};
    printf("%.3f\n", average_array(fill_array(array, SIZE_ARR), SIZE_ARR));

    return EXIT_SUCCESS;
}
#endif