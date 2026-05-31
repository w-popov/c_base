/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B21: Дан текст состоящий из английских букв и цифр, 
 *      заканчивается символом «.» 
 *      Перевести все заглавные английские буквы в строчные
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "HW5.h"

#ifndef TEST_DEF_HW5
int main (void)
{
    enum { STR_SIZE = 512 };
    char string[STR_SIZE] = {'\0'};
    char lower[STR_SIZE];
    scanf(" %[^.]", string);
    printf("%s\n", to_lowercase(string, lower, STR_SIZE));
    
    return EXIT_SUCCESS;
}
#endif


const char* to_lowercase (char* text, char* lower, const int str_size)
{
    memset(lower, '\0', str_size);
    enum { UPPER_SHIFT = 5 };
    
    uint8_t mask = 1 << UPPER_SHIFT;
    for (int i = 0; i < str_size; ++i)
    {
        if (!text[i])
            break;
        if (text[i] >= 0 && text[i] <= 32)
        {
            lower[i] = text[i];
            continue;
        }
        lower[i] = text[i] | mask;
    }
    return (const char*)lower;
}