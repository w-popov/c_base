/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D10: Дано натуральное число n ≥ 1. Проверьте, является ли оно простым.
 *      Программа должна вывести слово YES, если число простое или NO в
 *      противном случае. Необходимо составить рекурсивную функцию и
 * использовать ее.
 */
#include <stdio.h>
#include <stdlib.h>

int is_prime (int n, int delitel)
{
    if (n < 2)
    {
        return 0;
    }
    else if ((delitel * delitel) <= n)
    {
        if (!(n % delitel))
        {
            return 0;
        }
    }
    else
    {
        return 1;
    }

    return is_prime(n, delitel + 1);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%s\n", is_prime(input_number, 2) ? "YES" : "NO");
    return EXIT_SUCCESS;
}
#endif