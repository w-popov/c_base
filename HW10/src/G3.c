/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G3: В файле input.txt дана строка из не более 1000 символов. 
 *     Показать номера символов, совпадающих с последним символом строки.
 *     Результат записать в файл output.txt 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G3/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

const char* read_file_g3 (const char *file_name, char buffer[], const unsigned size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL; 
    }

    FILE *file = NULL;
    
    // "%1000[^\n]"
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

int last_character_number (const char buffer[], const unsigned int size, const char* file_name)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return 0; 
    }

    FILE *file = NULL;  
    // "%1000[^\n]"
    char format[12];
    snprintf(format, sizeof(format), "%%%d[^\n]", size - 1);
    
    size_t end_position = strlen(buffer) - 1;
    char last = buffer[end_position];
    
    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return 0;
    }

    int flag_space = 0;
    for (size_t i = 0; i < end_position; ++i)
    {
        if (last == buffer[i])
        {
            if (!flag_space)
            {
                fprintf(file, "%zu", i);
                flag_space = 1;
            }
            else 
            {
                fprintf(file, " %zu", i);
            }
        }
    }
    fprintf(file, "%s", "\n");
    fclose(file);
    return 1;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    if (read_file_g3(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE; 
    }
    if ( !last_character_number(buffer, MAX_LEN, FILE_OUTPUT) )
    {
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
#endif