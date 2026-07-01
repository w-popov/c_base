/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G17: В файле input.txt записаны символы. Необходимо разработать функцию,
 * которая меняет местами пары соседних символов не обращая внимание на символы
 * пробел. Если количество символов нечетно (пробелы не считаем), то последний
 * символ не меняем. Результат записать в файл output.txt.
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G17/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 1001

char *read_file_g17 (const char *file_name, char buffer[], const int size)
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

int write_to_file_g17 (const char *file_name, const char *string)
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

char *pairs_of_characters (char string[])
{
    if (string == NULL)
    {
        fprintf(stderr, "Пустой указатель\n");
        return NULL;
    }
    int len = strlen(string);
    int i = 0;
    /* Hello world! => eHllw orodl! */
    while (i < len - 1)
    {
        while (i < len && string[i] == ' ')
        {
            i++;
        }
        if (i >= len - 1)
        {
            break;
        }
        int first = i;
        i++;
        while (i < len && string[i] == ' ')
        {
            i++;
        }
        if (i >= len)
        {
            break;
        }
        int second = i;
        char temp = string[first];
        string[first] = string[second];
        string[second] = temp;
        i++;
    }
    return string;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char string[MAX_LEN] = {'\0'};
    if (read_file_g17(FILE_INPUT, string, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    char *s = pairs_of_characters(string);
    write_to_file_g17(FILE_OUTPUT, s);
    return EXIT_SUCCESS;
}
#endif