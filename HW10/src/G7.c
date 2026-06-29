/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G7: В файле input.txt считать символьную строку, не более 10 000 символов. 
 * Посчитать количество строчных (маленьких) и прописных (больших) букв в 
 * введенной строке. Учитывать только английские буквы. 
 * Результат записать в файл output.txt. 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G7/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 10001

struct UpperLower
{
    int upper;
    int lower;
};

char* read_file_g7 (const char *file_name, char buffer[], const int size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL; 
    }

    FILE *file = NULL;
    
    // "%10000[^\n]"
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

int write_to_file_g7 (const char *file_name, struct UpperLower ul)
{
    FILE *file = NULL;

    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return 0;
    }
    
    fprintf(file, "%d %d\n", ul.lower, ul.upper);
    fclose(file);
    return 1;
}

struct UpperLower number_of_letters (char buffer[])
{
    struct UpperLower ul = {.lower = 0, .upper = 0};
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return ul; 
    }
    for (size_t i = 0; i < strlen(buffer); ++i)
    {
        if (isalpha(buffer[i]))
        {
            if (isupper(buffer[i]))
            {
                ++(ul.upper);
            }
            else if (islower(buffer[i]))
            {
                ++(ul.lower);
            }
        }
    }
    return ul;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    
    if (read_file_g7(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    write_to_file_g7(FILE_OUTPUT, number_of_letters(buffer));
    return EXIT_SUCCESS;
}
#endif