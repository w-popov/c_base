/**
 * ДЗ-3. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * Ввести три числа, найти их среднее арифметическое.
 */

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    int a = 0, b = 0, c = 0;
    scanf("%d %d %d", &a, &b, &c);
    printf("%.2f\n", (float)(a + b + c) / 3); 

    return EXIT_SUCCESS;
}