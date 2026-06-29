/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G5: В файле input.txt дана символьная строка не более 1000 символов. 
 * Необходимо заменить все буквы "а" на буквы "b" и наоборот, как заглавные, 
 * так и строчные. Результат записать в output.txt. 
 */
#include <stdio.h>
#include <stdlib.h>

/*
  CATALOG_FILES_DEF_HW10 - полный путь к 
  каталогу файлов .txt у меня. 
  Определяется в HW10/src/CMakeLists.txt
*/
#ifdef CATALOG_FILES_DEF_HW10
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G5/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

char* read_file_g5 (const char *file_name, char buffer[], const int size)
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

void replace_A_with_B (char buffer[], const int size, const char *file_name)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return; 
    }
    for (int i = 0; i < size; ++i)
    {
        if (buffer[i] == 'a')
        {
            buffer[i] = 'b';
        }
        else if (buffer[i] == 'A')
        {
            buffer[i] = 'B';
        }
        else if (buffer[i] == 'b')
        {
            buffer[i] = 'a';
        }
        else if (buffer[i] == 'B')
        {
            buffer[i] = 'A';
        }
    }
    FILE *file = NULL;
    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return ;
    }
    fprintf(file, "%s\n", buffer);
    fclose(file);
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    if (read_file_g5(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    replace_A_with_B(buffer, MAX_LEN, FILE_OUTPUT);
    
    return EXIT_SUCCESS;
}
#endif