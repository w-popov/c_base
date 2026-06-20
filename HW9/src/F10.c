/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F10: Дана строка состоящая из маленьких латинских
 *      букв 'a'..'z'. В конце строки точка.
 *      Необходимо заменить повторяющиеся буквы на <буква><количество
 * повторений> Например: aaaaabbbc. ==> a5b3c1
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *pack_a_string (char *input_string, char *result, size_t result_size)
{
    if (scanf("%[^.]", input_string) != 1)
    {
        printf("Error scanf()\n");
        return input_string;
    }
    getchar();

    size_t len_string = strlen(input_string);
    if (len_string < 1)
    {
        return input_string;
    }

    char *destination = result;
    int count = 1;

    for (size_t i = 1; i <= len_string; ++i)
    {
        if (i < len_string && input_string[i] == input_string[i - 1])
        {
            count++;
        }
        else
        {
            size_t remaining = result_size - (destination - result);
            destination += snprintf(destination, remaining, "%c%d", input_string[i - 1], count);
            count = 1;
        }
    }
    *destination = '\0';

    return result;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { SIZE_ARR = 1000 };
    char array[SIZE_ARR] = {'\0'};
    char result[SIZE_ARR * 2] = {'\0'};
    printf("%s\n", pack_a_string(array, result, sizeof(result)));

    return EXIT_SUCCESS;
}
#endif