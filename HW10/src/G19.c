/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G19: Разработать функцию дана строка из маленьких английских букв.
 * Составить из символов палиндром максимальной длинны. При составлении
 * палиндрома буквы в палиндроме должны быть расположены в
 * лексикографическом порядке. Записать ответ в файл output.txt.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/*
  CATALOG_FILES_DEF_HW10 - полный путь к
  каталогу файлов .txt у меня.
  Определяется в HW10/src/CMakeLists.txt
*/
#ifdef CATALOG_FILES_DEF_HW10
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G19/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 1001

char *read_file_g19 (const char *file_name, char buffer[], const int size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL;
    }

    // "%1000[^\n]"
    char format[18];
    snprintf(format, sizeof(format), "%%%d[^\n]", size - 1);

    FILE *file = NULL;

    if ((file = fopen(file_name, "r")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return NULL;
    }
    fscanf(file, format, buffer);
    fclose(file);

    return buffer;
}

int write_to_file_g19 (const char *file_name, const char *string)
{
    FILE *file = NULL;

    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return 0;
    }
    fprintf(file, "%s\n", string);

    fclose(file);
    return 1;
}

void make_palindrome (char *string, char *result)
{
    enum { ALPHABET = 26 };
    int albet[ALPHABET] = {0};
    size_t len = strlen(string);

    /* счетчик символов */
    for (size_t i = 0; i < len; ++i)
    {
        if (string[i] >= 'a' && string[i] <= 'z')
        {
            albet[string[i] - 'a']++;
        }
    }

    size_t left_idx = 0;
    int center_char = -1;

    for (int i = 0; i < ALPHABET; ++i)
    {
        int pairs = albet[i] / 2;
        for (int j = 0; j < pairs; ++j)
        {
            result[left_idx++] = (char)('a' + i);
        }
        if (albet[i] % 2 != 0 && center_char == -1)
        {
            center_char = i;
        }
    }

    int len_left = (int)strlen(result) - 1;

    if (center_char != -1)
    {
        result[left_idx++] = (char)('a' + center_char);
    }
    
    /* копия реверса левой части вправо */
    for (int r = len_left; r >= 0; --r)
    {
        result[left_idx++] = result[r];
    }
    result[left_idx] = '\0';
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char string[MAX_LEN] = {'\0'};
    char palindrome[MAX_LEN] = {'\0'};
    read_file_g19(FILE_INPUT, string, MAX_LEN);
    make_palindrome(string, palindrome);
    write_to_file_g19(FILE_OUTPUT, palindrome);
    return EXIT_SUCCESS;
}
#endif