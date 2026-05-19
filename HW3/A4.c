/**
 * ДЗ-3. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * Ввести три числа, найти их сумму и произведение.
 */

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    int a = 0, b = 0, c = 0;
    scanf("%d %d %d", &a, &b, &c);
    printf(
        "%d+%d+%d=%d\n%d*%d*%d=%d\n", 
        a, b, c, (a + b + c), 
        a, b, c, (a * b * c)
    ); 

    return EXIT_SUCCESS;
}