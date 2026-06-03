/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C11: Составить функцию, которая определяет наибольший
 *      общий делитель двух натуральных и привести пример
 *      ее использования. int nod(int a, int b)
 */
#include <stdio.h>
#include <stdlib.h>

int nod (int a, int b)
{
    while (b)
    {
        int bucket = b;
        b = a % b;
        a = bucket;
    }
    return a > 0 ? a : -a;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_num_a = 0, input_num_b = 0;
    scanf("%d %d", &input_num_a, &input_num_b);
    printf("%d\n", nod(input_num_a, input_num_b));
    return EXIT_SUCCESS;
}
#endif