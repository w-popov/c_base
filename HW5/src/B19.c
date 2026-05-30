/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B8: Ввести натуральное число и определить, 
 *     верно ли, что сумма его цифр равна 10.
 */

#include <stdio.h>
#include <stdlib.h>

const char* is_summ_digits_ten(int number)
{
    enum { TEN = 10 };
    int sum = 0;
    if (number < 10)
        return "NO";
    while (number)
    {
        if (sum > TEN)
            break;
        sum += number % TEN;
        number /= TEN;
        if (sum == TEN)
            return "YES";
    }
    return "NO";    
}

#ifndef TEST_DEF_HW5
int main(void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%s\n", is_summ_digits_ten(number));

    return EXIT_SUCCESS;
}
#endif