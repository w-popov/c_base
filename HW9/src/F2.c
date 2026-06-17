/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F2: Написать только одну функцию, которая ставит в 
 *     начало массива все четные элементы, а в конец – все нечетные
 *     void sort_even_odd(int n, int a[]). 
 *     Не нарушайте порядок следования чисел между собой.
 */
#include <stdio.h>
#include <stdlib.h>

void sort_even_odd(int n, int a[])
{
    int *array = calloc(n, sizeof(int));
    int *array_even = array;
    
    for (int i = 0; i < n; ++i)
    {
        if (a[i] % 2 == 0)
        {
            *array_even++ = a[i];
        }
    }
    for (int i = 0; i < n; ++i)
    {
        if (a[i] % 2 != 0)
        {
            *array_even++ = a[i];
        }
    }
    for (int j = 0; j < n; ++j)
    {
        a[j] = array[j];
    }
    free(array);
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif