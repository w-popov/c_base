/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F20: Дан целочисленный массив из 10 элементов. 
 *      Необходимо определить количество четных и нечетных чисел.
 *      Если количество чётных чисел больше, чем количество нечётных,
 *      заменить каждое нечетное число на произведение нечетных цифр в 
 *      его десятичной записи. Если количество нечётных чисел больше или 
 *      равно количеству чётных, заменить каждое чётное число на 
 *      произведение чётных цифр в его десятичной записи.
 */
#include <stdio.h>
#include <stdlib.h>

enum { MORE_ODD, MORE_EVEN };

typedef int (*Compare)(int);

int is_even (int digit)
{
    return !(digit % 2);
}

int is_odd (int digit)
{
    return digit % 2;
}

/**
 * Произведение четных или нечетных цифр числа в зависимости от func.
 * @param func число четное или не четное?
 */
int mult_digits (int num, Compare func)
{
    if (!num)
    {
        return 1;
    }
    if (!func(num % 10))
    {
        return num % 10 * mult_digits(num / 10, func);
    }
    else
    {
        return mult_digits(num / 10, func);
    }
}

/**
 * Заполнить массив
 */
int *fill_array_F20 (int *array, const int size)
{
    if (array == NULL)
    {
        perror("Pointer array is NULL in func fill_array_F20()\n");
        return NULL;
    }
    for (int i = 0; i < size; ++i)
    {
        if (scanf("%d", &array[i]) != 1)
        {
            printf("Error scanf() in *fill_array_F20(int[], int)\n");
            return NULL;
        }
    }
    return array;
}

/**
 * Вернуть состояние в зависимости от того, каких (четных или нечетных)
 * чисел в массиве больше
 */
int get_more_even_or_odd (int *array, const int size)
{
    if (array == NULL)
    {
        perror("Pointer array is NULL in func get_even_nums()\n");
        return 0;
    }
    int evens = 0, odds = 0;
    for (int i = 0; i < size; ++i)
    {
        if (array[i] % 2)
        {
            odds++;
        }
        else
        {
            evens++;
        }
    }
    if (odds >= evens)
    {
        return MORE_ODD;
    }
    return MORE_EVEN;
}

/**
 * Произветси замену в массиве значений согласно условию,
 * в зависимости от четности/не четности.
 */
void operation (int *array, const int size, Compare func)
{
    if (array == NULL)
    {
        perror("Pointer array is NULL in func operation()\n");
        return;
    }
    for (int i = 0; i < size; ++i)
    {
        if (!func(array[i]))
        {
            if (array[i] == 0)
            {
                continue;
            }
            array[i] = mult_digits(array[i], func);
        }
    }
}

/**
 * Решение задания. (точка входа)
 */
int *even_and_odd (int *array, const int size)
{
    if (array == NULL)
    {
        perror("Pointer array is NULL in func even_and_odd()\n");
        return NULL;
    }
    int state = get_more_even_or_odd(array, size);
    switch (state)
    {
    case MORE_ODD:
        operation(array, size, is_odd);
        break;
    case MORE_EVEN:
        operation(array, size, is_even);
        break;
    default:
        break;
    }
    return array;
}

void print_array (int *array, const int size)
{
    if (array == NULL)
    {
        perror("Pointer array is NULL in func print_array()\n");
        return;
    }
    for (int i = 0; i < size; ++i)
    {
        printf("%d%s", array[i], i < size - 1 ? " " : "\n");
    }
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    print_array(
        even_and_odd(
            fill_array_F20(array, SIZE_ARR),
            SIZE_ARR
        ),
        SIZE_ARR
    );
    
    return EXIT_SUCCESS;
}
#endif