/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D8: Составить рекурсивную функцию, Выведите все
 *     числа от A до B включительно, в порядке возрастания,
 *     если A < B, или в порядке убывания в противном случае.
 */
#include <stdio.h>
#include <stdlib.h>

void print_from_A_to_B_rec (int a, int b)
{
    if (a < b)
    {
        printf("%d ", a);
        print_from_A_to_B_rec(a + 1, b);
    }
    else if (a > b)
    {
        printf("%d ", a);
        print_from_A_to_B_rec(a - 1, b);
    }
    else
    {
        printf("%d\n", a);
        return;
    }
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_a = 0, input_b = 0;
    scanf("%d %d", &input_a, &input_b);
    print_from_A_to_B_rec(input_a, input_b);
    return EXIT_SUCCESS;
}
#endif