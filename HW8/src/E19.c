/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E19: Вывести в порядке следования цифры, 
 *      входящие в десятичную запись натурального числа N 
 */
#include <stdio.h>
#include <stdlib.h>

void numbers_asc_order (int num)
{
    if (!num)
    {
        return;
    }
    numbers_asc_order(num / 10);
    printf("%d ", num % 10);
}

#ifndef TEST_DEF_HW8
int main (void)
{
    int input_number = 0;
    if (scanf("%d", &input_number) != 1)
    {
        printf("Error scanf\n");
        exit(EXIT_FAILURE);
    }
    numbers_asc_order(input_number);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif