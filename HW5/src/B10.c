/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B10: Ввести целое число и определить, верно ли, 
 *      что все его цифры расположены в порядке возрастания
 */

#include <stdio.h>
#include <stdlib.h>

const char* is_in_ascending_order(int number)
{
    // 123
    if (number < 10)
        return "YES";
    int prev_number = number % 10;
    do
    {
        number /= 10;
        if ((number % 10) >= prev_number)
            return "NO";
        prev_number = number % 10;
    } while (number);

    return "YES";
}


#ifndef TEST_DEF_HW5
int main(void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%s\n", is_in_ascending_order(number));

    return EXIT_SUCCESS;
}
#endif