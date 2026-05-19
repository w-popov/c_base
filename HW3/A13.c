/**
 * ДЗ-3. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * На вход подается произвольное трехзначное число,
 * напечать произведение цифр
 */

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    unsigned input_number = 0;
    unsigned x100 = 0, x10 = 0, x1 = 0;

    scanf("%u", &input_number);
    x100 = input_number / 100;
    x10 = (input_number % 100) / 10;
    x1 = input_number % 10;
    printf("%u\n", x100 * x10 * x1); 

    return EXIT_SUCCESS;
}