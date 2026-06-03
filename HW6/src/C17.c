/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C17: Составить логическую функцию, которая определяет,
 *      верно ли, что в заданном числе сумма цифр равна
 *      произведению. int is_happy_number(int n)
 */
#include <stdio.h>
#include <stdlib.h>

int is_happy_number (int n)
{
    enum { ZERO, ONE, TEN = 10 };
    int sum = ZERO, mul = ONE;
    while (n)
    {
        int rest = n % TEN;
        sum += rest;
        mul *= rest;
        n /= TEN;
    }
    return sum == mul ? ONE : ZERO;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int input_num = 0;
    scanf("%d", &input_num);
    printf("%s\n", is_happy_number(input_num) ? "YES" : "NO");
    return EXIT_SUCCESS;
}
#endif