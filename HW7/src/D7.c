/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D7: Составить рекурсивную функцию, печать всех чисел от N до 1.
 */
#include <stdio.h>
#include <stdlib.h>

void print_from_N_to_1_rec (int num)
{
    if (!num)
    {
        return;
    }
    printf("%d ", num);
    print_from_N_to_1_rec(num - 1);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    print_from_N_to_1_rec(input_number);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif