/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F8: В последовательности записаны целые числа от M до N 
 *     ( M меньше N, M больше или равно 1)в произвольном порядке, 
 *     но одно из чисел пропущено (остальные встречаются ровно по одному разу).
 *     N не превосходит 1000. Последовательность заканчивается числом 0. 
 *     Определить пропущенное число. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int missing_number (int *array, const int size)
{
    int *cp_array = malloc(size * sizeof(int));
    if (cp_array == NULL)
    {
        printf("Error alloc mem in missing_number(int[], const int)\n");
        return -1;
    }
    /* Копия массива */
    memcpy(cp_array, array, size * sizeof(int));

    /* Сортировка вставками */
    for (int i = 1; i < size; ++i)
    {
        int key = cp_array[i];
        int j = i - 1;
        while (j >= 0 && cp_array[j] > key)
        {
            cp_array[j + 1] = cp_array[j];
            j--;
        }
        cp_array[j + 1] = key;
    }
    /* Поиск */
    for (int i = 1; i < size; ++i)
    {
        if (cp_array[i] != (cp_array[i-1] + 1))
        {
            return cp_array[i] - 1;
        }
    }
    
    free(cp_array);
    return 0;
}

int fill_array_F8 (int *array, const int size)
{
    int size_full = 0;
    for (int i = 0; i < size; ++i)
    {
        int number = 0;
        int return_scanf = scanf("%d", &number);
        if (return_scanf < 1)
        {
            return 0;
        }
        size_full += return_scanf;
        if (!number)
        {
            --size_full;
            break;
        }
        array[i] = number;
    }
    return size_full;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { SIZE_ARR = 1000 };
    int array[SIZE_ARR] = {0};
    const int full_size = fill_array_F8(array, SIZE_ARR);
    printf("%d\n", missing_number(array, full_size));
    return EXIT_SUCCESS;
}
#endif