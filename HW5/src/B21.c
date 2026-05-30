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

enum Sizies { UPPER_SHIFT=5, STR_SIZE=512 };

const char* to_lowercase (char* text, char* lower)
{
    memset(lower, '\0', STR_SIZE);
    
    uint8_t mask = 1 << UPPER_SHIFT;
    for (int i = 0; i < STR_SIZE; ++i)
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

#ifndef TEST_DEF_HW5
int main(void)
{
    char string[STR_SIZE] = {'\0'};
    char lower[STR_SIZE];
    scanf(" %[^.]", string);
    printf("%s\n", to_lowercase(string, lower));

    return EXIT_SUCCESS;
}
#endif