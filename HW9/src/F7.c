/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F7: Написать функцию, которая сжимает серии массива,
 * состоящего из единиц и нулей по следующему принципу: например,
 * массив [0,0,0,0,1,1,1,1,1,1,1,0,0,1,1,1,1] => [4,7,2,4]
 * а массив [1,1,1,0,0,0,0,0,0,0] преобразуется в [0,3,7]
 * Прототип функции: int compression(int a[], int b[], int N)
 */
#include <stdio.h>
#include <stdlib.h>

/**
 * Функция принимает исходный массив a[]
 * и сжимает в массив b[], N - число элементов в массиве a[].
 */
int compression (int a[], int b[], int N)
{
    if (!a || !b || N <= 0)
    {
        return 0;
    }

    int count_b = 0;
    if (a[0] == 1)
    {
        b[count_b++] = 0;
    }

    int current_len = 1;

    for (int i = 1; i < N; ++i)
    {
        if (a[i] == a[i - 1])
        {
            ++current_len;
        }
        else
        {
            b[count_b++] = current_len;
            current_len = 1;
        }
    }
    
    b[count_b++] = current_len;

    return count_b;
}

#ifndef TEST_DEF_HW9
int main (void)
{

    return EXIT_SUCCESS;
}
#endif