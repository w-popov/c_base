/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A14: Дано трехзначное положительное целое число,
 * напечатать макисмальную цифру [123 => 3]
 */

#include <stdio.h>
#include <stdlib.h>

/**
 * Вернуть макисмальную цифру [123 => 3]
 */
unsigned max_digit_from_number(unsigned input_number)
{
    unsigned x100 = 0, x10 = 0, x1 = 0;
    x100 = input_number / 100;
    x10 = (input_number % 100) / 10;
    x1 = input_number % 10;
    
    return (x100 > x10) ? ((x100 > x1) ? x100 : x1) : ((x10 > x1) ? x10 : x1);
}

/**
 * TEST_DEF_HW определяется в Makefile для компиляции тестов в ./tests. 
 * Компиляция без main()
 */
#ifndef TEST_DEF_HW
int main(void)
{
    unsigned input_number = 0;
    scanf("%u", &input_number);
    printf("%u\n", max_digit_from_number(input_number));

    return EXIT_SUCCESS;
}
#endif