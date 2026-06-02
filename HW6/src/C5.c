/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C5: Составить функцию, которая определяет сумму
 *     всех чисел от 1 до N и привести пример ее использования
 */
#include <stdio.h>
#include <stdlib.h>

int summ_from_one_to_number (int target_number)
{
    int summ = 0;
    for (int i = 1; i <= target_number; ++i)
    {
        summ += i;
    }
    return summ;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%d\n", summ_from_one_to_number(number));
    return EXIT_SUCCESS;
}
#endif