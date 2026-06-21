/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F11: Необходимо создать функцию, которая находит и
 *      выводит в порядке возрастания номера двух элементов
 *      массива, сумма которых минимальна
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

struct MinTwo
{
    int min1;
    int min2;
};

int *fill_array_F11 (int *array, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        if (scanf("%d", &array[i]) != 1)
        {
            printf("Error scanf() in *fill_array_F11(int[], int)\n");
            return NULL;
        }
    }
    return array;
}

struct MinTwo find_2_minimum (int *array, const int size)
{
    if (array == NULL)
    {
        printf("Error array pointer in struct MinTwo find_2_minimum(int[], int)\n");
        return (struct MinTwo){.min1 = -1, .min2 = -1};
    }
    // 1 2
    int index_min1 = 0;
    for (int i = 1; i < size; ++i)
    {
        if (array[i] < array[index_min1])
        {
            index_min1 = i;
        }
    }

    int index_min2 = (!index_min1) ? 1 : 0;
    for (int i = 0; i < size; ++i)
    {
        if ((i != index_min1) && (array[i] < array[index_min2]))
        {
            index_min2 = i;
        }
    }

    struct MinTwo mt;
    
    mt.min1 = index_min1 < index_min2 ? index_min1 : index_min2;
    mt.min2 = index_min2 > index_min1 ? index_min2 : index_min1;
    return mt;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { SIZE_ARR = 30 };
    int array[SIZE_ARR] = {0};
    if (fill_array_F11(array, SIZE_ARR) == NULL)
    {
        return EXIT_FAILURE;
    }
    struct MinTwo mt = find_2_minimum(array, SIZE_ARR);
    printf("%d %d\n", mt.min1, mt.min2);
    return EXIT_SUCCESS;
}
#endif