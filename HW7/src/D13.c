/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D13: Составить рекурсивную функцию, печать всех простых множителей числа
 */
#include <stdio.h>
#include <stdlib.h>

void print_simple_rec (int num, int divider)
{
    if (num < 2)
    {
        return;
    }
    if (!(num % divider))
    {
        printf("%d ", divider);
        print_simple_rec(num / divider, divider);
    }
    else
    {
        if (divider * divider > num)
        {
            printf("%d", num);
            return;
        }
        print_simple_rec(num, divider + 1);
    }
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    print_simple_rec(input_number, 2);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif