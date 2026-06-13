/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E15: Считать 10 чисел в диапазоне от -500 до 500 и 
 *      разложить по двум массивам: в одни помещать только 
 *      положительные, во второй - только отрицательные. 
 *      Числа, равные нулю, игнорировать. 
 *      Вывести на экран все элементы обоих массивов. 
 */
#include <stdio.h>
#include <stdlib.h>

void print_array (int *pos_arr, int *neg_arr, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        if (pos_arr[i])
        {
            printf("%d ", pos_arr[i]);
        }
    }
    for (int i = 0; i < size_array; ++i)
    {
        if (neg_arr[i])
        {
            printf("%d ", neg_arr[i]);
        }
    }
    printf("\n");
}

void positive_negative_arrs (
                                int *array,
                                int *pos_arr,
                                int *neg_arr,
                                const int size
                            )
{
    for (int i = 0, p = 0, n = 0; i < size; ++i)
    {
        if (array[i])
        {
            if (array[i] > 0)
            {
                pos_arr[p++] = array[i];
            }
            else
            {
                neg_arr[n++] = array[i];
            }
        }
    }
}

int *fill_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        scanf("%d", &array[i]);
    }
    return array;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    int neg_array[SIZE_ARR] = {0};
    int pos_array[SIZE_ARR] = {0};

    positive_negative_arrs(
        fill_array(array, SIZE_ARR), 
        pos_array,
        neg_array,
        SIZE_ARR
    );
    print_array(pos_array, neg_array, SIZE_ARR);    
    
    return EXIT_SUCCESS;
}
#endif