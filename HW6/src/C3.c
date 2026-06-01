/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C3: Написать функцию, которая возвращает среднее
 *     арифметическое двух переданных ей аргументов
 *     (параметров). int middle(int a, int b)
 */
#include <stdio.h>
#include <stdlib.h>

int middle (int, int);

#ifndef TEST_DEF_HW6
int main (void)
{
    int number_1 = 0, number_2 = 0;
    scanf("%d %d", &number_1, &number_2);
    printf("%d\n", middle(number_1, number_2));
    return EXIT_SUCCESS;
}
#endif

int middle (int number1, int number2) { return (number1 + number2) / 2; }