/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C15: Составить функцию логическую функцию, которая определяет,
 * верно ли, что в заданном числе все цифры стоят по возрастанию.
 * Используя данную функцию решить задачу. int grow_up(int n)
 */
#include <stdio.h>
#include <stdlib.h>

int grow_up (int n)
{
    enum { ZERO, ONE, TEN = 10 };
    if (n < TEN)
    {
        return ONE;
    }
    while (n)
    {
        int rest = n % TEN;
        n /= TEN;
        if (n % TEN >= rest)
        {
            return ZERO;
        }
    }
    return ONE;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_num = 0;
    scanf("%d", &input_num);
    printf("%s\n", grow_up(input_num) ? "YES" : "NO");
    return EXIT_SUCCESS;
}
#endif