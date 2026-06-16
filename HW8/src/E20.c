/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E20: Переставить цифры в числе так, чтобы 
 *      получилось максимальное число 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef int (*Compare)(int, int);

int compare_desc (int a, int b)
{
    return a > b;
}

/**
 * Сортировка массива int[] вставками
 */
int *sort_X_array (int *array, const int size_array, Compare cmp)
{
    for (int i = 1; i < size_array; ++i)
    {
        int key = array[i];
        int j = i - 1;
        while (j >= 0 && cmp(key, array[j]))
        {
            array[j + 1] = array[j];
            j--;
        }
        array[j + 1] = key;
    }
    return array;
}

int build_greatest_number (int *digits_array, int size_num)
{
    int *sort_array = sort_X_array(digits_array, (const int)size_num, compare_desc);
    int greatest = 0;
    for (int i = 0; i < size_num; ++i)
    {
        greatest = greatest * 10 + sort_array[i];  
    }
    return greatest;
}

int the_greatest_number (int number)
{
    if (!number)
    {
        return 0;
    }
    int greatest_number = 0, size_num = 0;
    int copy_number = number;
    while (copy_number)
    {
        copy_number /= 10;
        ++size_num;
    }
    int *array = calloc(size_num, sizeof(int));
    int *digits_array = array;
    while (number)
    {
        *(digits_array++) = number % 10;
        number /= 10;
    }
    greatest_number = build_greatest_number(array, size_num);
    /* free */
    free(array);

    return greatest_number;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    int input_number = 0;
    if (scanf("%d", &input_number) != 1)
    {
        printf("Error scanf\n");
        exit(EXIT_FAILURE);
    }
    printf("%d\n", the_greatest_number(input_number));
    return EXIT_SUCCESS;
}
#endif