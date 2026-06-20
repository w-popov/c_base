/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F6: Написать только одну логическую функцию, которая определяет, верно ли, 
 *     что среди элементов массива есть два одинаковых. Если ответ «да», функция возвращает 1; 
 *     если ответ «нет», то 0. Строго согласно прототипу:
 *     int is_two_same(int size, int a[]);
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int is_two_same(int size, int a[])
{
    int *copy_a = malloc(size * sizeof(int));
    if (copy_a == NULL)
    {
        printf("Error alloc mem in is_two_same(int, int[])\n");
        return -1;
    }
    /* Копия массива */
    memcpy(copy_a, a, size * sizeof(int));

    /* Сортировка вставками */
    for (int i = 1; i < size; ++i)
    {
        int key = copy_a[i];
        int j = i - 1;
        while (j >= 0 && copy_a[j] > key)
        {
            copy_a[j + 1] = copy_a[j];
            j--;
        }
        copy_a[j + 1] = key;
    }
    int is_identical = 0;
    for (int i = 1; i < size; ++i)
    {
        if (copy_a[i] == copy_a[i - 1])
        {
            is_identical = 1;
            break;
        }
    }
    free(copy_a);
    return is_identical;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    return EXIT_SUCCESS;
}
#endif