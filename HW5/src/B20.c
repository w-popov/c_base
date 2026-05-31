/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B20: Проверить число на простоту.
 */

#include <stdio.h>
#include <stdlib.h>

const char* is_prime_number (int number)
{
    if (number < 2)
        return "NO";
    for (int counter = 2; (counter * counter) <= number; ++counter) {
        if ( !(number % counter) ) 
            return "NO";
    }
    return "YES";
}

#ifndef TEST_DEF_HW5
int main (void)
{
    int number = 0;
    scanf("%d", &number);
    printf("%s\n", is_prime_number(number));

    return EXIT_SUCCESS;
}
#endif