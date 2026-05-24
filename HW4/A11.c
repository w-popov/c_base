/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A11: Напечатать сумму максимума и минимума.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef _Bool(*compare)(int, int);

_Bool x_great_then_y(int x, int y)
{
    return x > y;
}

_Bool x_less_then_y(int x, int y)
{
    return x < y;
}

void compare_nums
                (
                    int b, 
                    int c, 
                    int d, 
                    int e, 
                    int* target, 
                    compare compare_func
                )
{
    if (!target)
        exit(EXIT_FAILURE);
    if (compare_func(b, *target)) *target = b;
    if (compare_func(c, *target)) *target = c;
    if (compare_func(d, *target)) *target = d;
    if (compare_func(e, *target)) *target = e;

}

/**
 * Вернуть сумму максимума и минимума из 5 чисел.
 */
int sum_min_max_from_5_numbers(int a, int b, int c, int d, int e)
{
    int max = a, min = a;
    compare_nums(b, c, d, e, &max, x_great_then_y);
    compare_nums(b, c, d, e, &min, x_less_then_y);
        
    return min + max;
}


/**
 * TEST_DEF_HW определяется в Makefile для компиляции тестов в ./tests. 
 * Компиляция без main()
 */
#ifndef TEST_DEF_HW
int main(void)
{
    int a = 0, b = 0, c = 0, d = 0, e = 0;
    scanf("%d %d %d %d %d", &a, &b, &c, &d, &e);
    printf("%d", sum_min_max_from_5_numbers(a, b, c, d, e));

    return EXIT_SUCCESS;
}
#endif