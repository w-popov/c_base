/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C4: Описать функцию вычисления f(x) по формуле:
 * f(x)= x*x при -2 ≤ x < 2;
 * x*x + 4x + 5 при x ≥ 2;
 * 4 при x < -2.
 * Используя эту функцию для n заданных чисел, вычислить f(x).
 * Среди вычисленных значений найти наибольшее
 */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

typedef int (*function)(int);

/**
 * f(x)= x*x при -2 ≤ x < 2
 */
int square_function (int arg)
{
    return arg * arg;
}

/**
 * f(x) = x*x + 4x + 5 при x ≥ 2
 */
int square_solver_function (int arg)
{
    return (arg * arg) + (4 * arg) + 5;
}

/**
 * f(x) = 4 при x < -2
 */
int get_four_function (int arg)
{
    enum { FOUR = 4 };
    return FOUR;
}

/**
 * Вычислить функцию по аргументу
 */
int calculate_argument (int arg)
{
    function func = NULL;
    if (arg < -2)
    {
        func = get_four_function;
    }
    else if ((arg >= -2) && (arg < 2))
    {
        func = square_function;
    }
    else
    {
        func = square_solver_function;
    }

    if (func != NULL)
    {
        return func(arg);
    }
    return INT32_MAX;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    enum { SIZE_STR = 512 };
    int accumulate = INT32_MAX;
    int max = INT32_MIN;
    while (accumulate && scanf("%d", &accumulate))
    {
        int calculate = calculate_argument(accumulate);
        max = calculate > max ? calculate : max; 
    }
    printf("%d\n", max);
    return EXIT_SUCCESS;
}
#endif