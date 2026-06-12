/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E6: Считать массив из 12 элементов и посчитать 
 *     среднее арифметическое элементов массива 
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

float average_2_array (int *array, const int size_array)
{
    float avg = 0.0f;
    for (int i = 0; i < size_array; ++i)
    {
        avg += (float)(array[i]);
    }
    return avg / size_array;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 12 };
    int array[SIZE_ARR] = {0};
    printf("%.2f\n",
           average_2_array(
            fill_array(array, SIZE_ARR), 
            SIZE_ARR)
        );
    
    return EXIT_SUCCESS;
}
#endif