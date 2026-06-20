/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F9: Составить только функцию которая в массиве находит максимальный из 
 *     отрицательных элементов и меняет его местами с последним элементом массива. 
 *     Гарантируется, что в массиве только один такой элемент или же такой элемент 
 *     отсутствует. Если отрицательных элементов нет - массив не менять.
 *     Прототип функции: void swap_negmax_last(int size, int a[])
 */
#include <stdio.h>
#include <stdlib.h>

void swap_negmax_last(int size, int a[])
{
    #define INT32_MIN (-2147483647-1)
    if (size < 2)
    {
        return;
    }
    int last_index = size - 1;
    int max_min = INT32_MIN, index_max_min = 0;
    for (int i = 0; i < size; ++i)
    {
        if ((a[i] < 0) && (a[i] > max_min))
        {
            max_min = a[i];
            index_max_min = i;  
        }
    }
    if (max_min != INT32_MIN)
    {
        int temp_value = a[index_max_min];
        a[index_max_min] = a[last_index]; 
        a[last_index] = temp_value;
    }
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif