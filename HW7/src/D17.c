/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D17: Реализуйте функцию Аккермана по прототипу
 *      int akkerman(int m, int n)
 */
#include <stdio.h>
#include <stdlib.h>

int akkerman (int m, int n)
{
    if (!m)
    {
        return n + 1;
    }
    if (m && !n)
    {
        return akkerman(m - 1, 1);
    }
    return akkerman(m - 1, akkerman(m, n - 1));
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number_m = 0, input_number_n = 0;
    scanf("%d %d", &input_number_m, &input_number_n);
    printf("%d\n", akkerman(input_number_m, input_number_n));
    return EXIT_SUCCESS;
}
#endif