/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F3: Вывести в порядке возрастания цифры, входящие в
 *     десятичную запись натурального числа N, не более 1000 цифр.
 *     Цифра пробел сколько раз данная цифра встречается в числе.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

struct Digit_t
{
    uint16_t digit : 4;
    uint16_t count : 12;
};

/**
 * ВЫДЕЛЯЕТ ПАМЯТЬ!
 * Количество цифр в числе в виде массива структур 
 * 0: [0..кол-во цифр]
 * 1: [0..кол-во цифр]
 * . . .
 * 9: [0..кол-во цифр]
 */
struct Digit_t* numbers_in_asc_order (const char *digits, size_t num_size, const int max10)
{
    struct Digit_t *counter = calloc(max10, sizeof(struct Digit_t));
    for (int d = 0; d < max10; ++d)
    {
        counter[d].digit = d;
        counter[d].count = 0;
    }
    for (size_t i = 0; i < num_size; ++i)
    {
        char c = digits[i];
        if (c >= '0' && c <= '9')
        {
            counter[c - '0'].count++;
        }
    }
    return counter;
}

void printF3 (struct Digit_t* counter, const int max10)
{
    for (int i = 0; i < max10; ++i)
    {
        if (counter[i].count)
        {
            printf("%d %d", counter[i].digit, counter[i].count);
            printf("%s", i < max10 ? "\n" : "");
        }
    }
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { MAX10 = 10, SIZE = 1001 };
    char input_num[SIZE] = {0};
    if (scanf("%s", &input_num[0]) != 1)
    {
        return EXIT_FAILURE;
    }
    struct Digit_t* counter = numbers_in_asc_order(
                                input_num, 
                                strlen(input_num),
                                MAX10
                            );
    printF3(counter, MAX10);
    free(counter);
    return EXIT_SUCCESS;
}
#endif