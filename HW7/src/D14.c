/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D14: Дана последовательность целых чисел через пробел,
 *      завершающаяся числом 0. Выведите все НЕЧЕТНЫЕ числа
 *      ИЗ этой последовательности, сохраняя их порядок
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void print_odd_numbers_rec (int *pt_eq)
{
    if (pt_eq == NULL || !(*pt_eq))
    {
        return;
    }
    if (*pt_eq % 2)
    {
        printf("%d ", *pt_eq);
    }
    print_odd_numbers_rec(pt_eq + 1);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    enum { SIZE_EQ = 1024 };
    int input_eq[SIZE_EQ] = {0};
    int input_num = INT32_MIN;
    for (int i = 0; input_num && i < SIZE_EQ; ++i)
    {
        scanf("%d", &input_num);
        input_eq[i] = input_num;
    }
    print_odd_numbers_rec(&input_eq[0]);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif