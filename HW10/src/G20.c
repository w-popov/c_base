/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G20: Считать предложение из файла input.txt и определить можно ли из
 * английских букв предложения записанного в файле получить
 * одно слово - палиндром. Ответ напечатать на стандартный поток вывода.
 * Требуется реализовать логическую функцию и применить ее. is_palindrom(string)
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G20/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 1001

char *read_file_g20 (const char *file_name, char buffer[], const int size)
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

/*
 * палиндром возможен, если нечетных букв 0 или 1
 */
int is_palindrom (char string[])
{
    enum { ALPHABET = 26 };
    int albet[ALPHABET] = {0};
    int has_letters = 0;
    size_t len = strlen(string);

    /* счетчик символов */
    for (size_t i = 0; i < len; ++i)
    {
        if (string[i] >= 'a' && string[i] <= 'z')
        {
            albet[string[i] - 'a']++;
            has_letters = 1;
        }
    }
    if (!has_letters)
    {
        return 0;
    }
    int odd_count = 0;
    for (int i = 0; i < ALPHABET; ++i)
    {
        if (albet[i] % 2 != 0)
        {
            odd_count++;
        }
    }
    return (int)(odd_count <= 1);
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char string[MAX_LEN] = {'\0'};
    read_file_g20(FILE_INPUT, string, MAX_LEN);
    printf("%s\n", is_palindrom(string) ? "YES" : "NO");
    return EXIT_SUCCESS;
}
#endif