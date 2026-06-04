/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D4: Дано натуральное число N. Выведите все его цифры
 *      по одной, в прямом порядке, разделяя их пробелами
 *      или новыми строками. Необходимо реализовать рекурсивную функцию.
 *      void print_num(int num)
 */
#include <stdio.h>
#include <stdlib.h>
// #include "HW7.h"

void print_num (int num)
{
    if (num / 10)
    {
        print_num(num / 10);
    }
    printf("%d ", num % 10);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    print_num(input_number);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif