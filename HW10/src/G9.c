/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G9: В файле input.txt строка из меленьких и больших английских букв, 
 * знаков препинания и пробелов. Требуется удалить из нее повторяющиеся символы
 *  и все пробелы. Результат записать в файл output.txt. 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G9/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

char* read_file_g9 (const char *file_name, char buffer[], const int size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL; 
    }

    FILE *file = NULL;
    
    // "%1000[^\n]"
    char format[18];
    snprintf(format, sizeof(format), "%%%d[^\n]", size - 1);

    if ((file = fopen(file_name, "r")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return NULL;
    }
    if (fscanf(file, format, buffer) < 0)
    {
        fprintf(stderr, "Ошибка чтения или пустой файл: %s\n", file_name);
        fclose(file);
        return NULL;
    }
    fclose(file);
    
    return buffer;
}

int write_to_file_g9 (const char *file_name, const char string[])
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

char *remove_duplicate_chars (char buffer[], char str[])
{
    if (buffer == NULL || str == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[] или str[]\n");
        return NULL;
    }
    char *ptstr = str;
    int is_duplicate = 0;
    for (size_t i = 0; i < strlen(buffer); ++i)
    {
        is_duplicate = 0;
        if(isspace(buffer[i]))
        {
            continue;
        }
        for (size_t j = 0; j < strlen(buffer); ++j)
        {
            if (buffer[i] == str[j])
            {
                is_duplicate = 1;
                break;
            }
        }
        if (!is_duplicate)
        {
            *ptstr++ = buffer[i];
        }
    }
    return str;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    char cleanstr[MAX_LEN] = {'\0'};
    
    if (read_file_g9(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    if (remove_duplicate_chars(buffer, cleanstr) == NULL)
    {
        return EXIT_FAILURE;
    }
    write_to_file_g9(FILE_OUTPUT, cleanstr);
    
    return EXIT_SUCCESS;
}
#endif