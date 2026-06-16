/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E18: В диапазоне натуральных чисел от 2 до N определить, 
 *      сколько из них кратны любому из чисел в диапазоне от 2 до 9 
 */
#include <stdio.h>
#include <stdlib.h>

void multiplicity_from_2_to_99 (int num)
{
    int quants = 0;
    for (int q = 2; q <= 9; ++q)
    {
        for (int n = 2; n <= num; ++n)
        {
            if (!(n % q))
            {
                ++quants;
            }
        }
        printf("%d %d\n", q, quants);
        quants = 0;
    }
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
    multiplicity_from_2_to_99(input_number);
    return EXIT_SUCCESS;
}
#endif