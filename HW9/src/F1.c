/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F1: Написать только одну функцию, которая сортирует массив по возрастанию
 *     void sort_array(int size, int a[])
 */
#include <stdio.h>
#include <stdlib.h>

/* Сортировка вставками */
void sort_array(int size, int a[])
{
    for (int i = 1; i < size; ++i)
    {
        int key = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > key)
        {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = key;
    }
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { SIZE_ARR = 10 };
    // int array[SIZE_ARR] = {0};
    
    return EXIT_SUCCESS;
}
#endif