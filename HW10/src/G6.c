/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G6: В файле input.txt символьная строка не более 1000 символов. 
 * Необходимо проверить, является ли она палиндромом 
 * (палиндром читается одинаково в обоих направлениях). 
 * Реализовать логическую функцию is_palindrom(str) и 
 * записать ответ в файл output.txt.
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G6/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

char* read_file_g6 (const char *file_name, char buffer[], const int size)
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

int write_to_file_g6 (const char *file_name, const char *buffer)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return 0; 
    }

    FILE *file = NULL;

    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return 0;
    }
    
    fprintf(file, "%s\n", buffer);
    fclose(file);
    return 1;
}

int is_palindrom (char buffer[], const int size)
{
    int is_string_palindrom = 0;
    char *palindrom = calloc(size, sizeof(char));
    if (palindrom == NULL || buffer == NULL)
    {
        return 0;
    }
    size_t stlen = strlen(buffer);
    for (size_t i_pal = stlen - 1, i_buf = 0; i_buf < stlen; --i_pal, ++i_buf)
    {
        palindrom[i_pal] = buffer[i_buf];
    }
    if (!strncmp(buffer, palindrom, stlen))
    {
        is_string_palindrom = 1;
    }
    free(palindrom);
    return is_string_palindrom;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    const char *yes = "YES", *no = "NO";
    if (read_file_g6(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    if (is_palindrom(buffer, MAX_LEN))
    {
        write_to_file_g6(FILE_OUTPUT, yes);
    }
    else
    {
        write_to_file_g6(FILE_OUTPUT, no);
    }
    return EXIT_SUCCESS;
}
#endif