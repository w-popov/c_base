/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D16: Написать логическую рекурсивную функцию и используя ее
 *      определить является ли введенное натуральное число точной степенью
 *      двойки. int is2pow(int n)
 */
#include <stdio.h>
#include <stdlib.h>

int is2pow (int n)
{
    if (n <= 0)
    {
        return 0;
    }
    if (n == 1)
    {
        return 1;
    }
    if (n % 2)
    {
        return 0;
    }

    return is2pow(n / 2);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%s\n", is2pow(input_number) ? "YES" : "NO");
    return EXIT_SUCCESS;
}
#endif