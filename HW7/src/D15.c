/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D15: Дана последовательность ненулевых целых чисел,
 *      завершающаяся числом 0. Ноль в последовательность не входит.
 *      Определите наибольшее значение числа в этой последовательности.
 *      Для решения задачи необходимо написать рекурсивную функцию вида:
 *      int max_find(int max)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int max_find (int max)
{
    int num = INT32_MIN;
    if (scanf("%d", &num) != 1 || !num)
    {
        return max;
    }
    if (num > max)
    {
        max = num;
    }
    return max_find(max);
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%d\n", max_find(input_number));
    return EXIT_SUCCESS;
}
#endif