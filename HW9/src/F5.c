/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F5: Написать только одну функцию, которая находит 
 *     максимальный элемент в массиве. Всю программу загружать не надо.
 *     Прототип функции: int find_max_array(int size, int a[])
 */
#include <stdio.h>
#include <stdlib.h>

int find_max_array(int size, int a[])
{
    #define INT32_MIN (-2147483647-1)
    int max = INT32_MIN;
    for (int i = 0; i < size; ++i)
    {
        max = a[i] > max ? a[i] : max;
    }
    return max;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif