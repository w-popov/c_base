/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G1: В файле input.txt дана строка. Вывести ее в 
 *     файл output.txt три раза через запятую и 
 *     показать количество символов в ней.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
  CATALOG_FILES_DEF_HW10 - полный путь к 
  каталогу файлов .txt у меня. 
  Определяется в HW10/src/CMakeLists.txt
*/
#ifdef CATALOG_FILES_DEF_HW10
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G1/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 101

char* read_file (const char *file_name, char buffer[], const unsigned size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL; 
    }

    FILE *file = NULL;
    
    // "%99[^\n]"
    char format[12];
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

int write_to_file (const char *file_name, char buffer[], const unsigned size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return 0; 
    }

    FILE *file = NULL;
    enum { REPEAT = 3 };
    
    // "%99[^\n]"
    char format[12];
    snprintf(format, sizeof(format), "%%%d[^\n]", size - 1);
    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return 0;
    }
    for (int i = 0; i < REPEAT; ++i)
    {
        fprintf(file, "%s%s", buffer, i < REPEAT - 1 ? ", " : " ");
    }
    fprintf(file, "%zu\n", strlen(buffer));
    fclose(file);
    return 1;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    if (read_file(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        perror("Error in main\n");
        return EXIT_FAILURE;
    }
    write_to_file(FILE_OUTPUT, buffer, MAX_LEN);
    return EXIT_SUCCESS;
}
#endif