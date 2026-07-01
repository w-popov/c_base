/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G18: В файле input.txt необходимо удалить все лишние 
 * пробелы (в начале предложения и сдвоенные пробелы). 
 * Для решения задачи разработать функцию. 
 * Результат записать в output.txt. 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G18/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 1001

char *read_file_g18 (const char *file_name, char buffer[], const int size)
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

int write_to_file_g18 (const char *file_name, const char *string)
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

char *remove_spaces(char *string)
{
    if (string == NULL)
    {
        fprintf(stderr, "Пустой указатель\n");
        return NULL;
    }

    char *read = string;
    char *write = string;

    while (isspace((unsigned char)*read)) 
    {
        read++;
    }

    while (*read != '\0') 
    {
        while (*read != '\0' && !isspace((unsigned char)*read)) 
        {
            *write++ = *read++;
        }
        int has_space = 0;
        while (isspace((unsigned char)*read)) 
        {
            has_space = 1;
            read++;
        }
        if (has_space && *read != '\0') 
        {
            *write++ = ' ';
        }
    }
    *write = '\0'; 

    return string;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char string[MAX_LEN] = {'\0'};
    if (read_file_g18(FILE_INPUT, string, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    const char *s = remove_spaces(string);
    write_to_file_g18(FILE_OUTPUT, s);    
    return EXIT_SUCCESS;
}
#endif