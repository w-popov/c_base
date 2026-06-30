/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G10: В файле input.txt дана строка слов, разделенных пробелами.
 * Найти самое длинное слово и вывести его в файл output.txt. 
 * Случай, когда самых длинных слов может быть несколько, не обрабатывать 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/*
  CATALOG_FILES_DEF_HW10 - полный путь к 
  каталогу файлов .txt у меня. 
  Определяется в HW10/src/CMakeLists.txt
*/
#ifdef CATALOG_FILES_DEF_HW10
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G10/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

char* read_file_g10 (const char *file_name, char buffer[], const int size)
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

int write_to_file_g10 (const char *file_name, const char string[])
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

/**
 * @attention ВЫДЕЛЯЕТ дин. память!
 * @brief Разбить строку на части по разделителю и вернуть слово максимальной длины.
 * @details strspn - Возвращает длину начального сегмента str, 
 *                   состоящего только из символов, перечисленных в accept.
 *          strcspn - Возвращает длину начального сегмента str, 
 *                    не содержащего ни одного символа из reject
 * @param string входная строка.
 * @param delimeter разделитель.
 * @return Строка максимальной длины.
 * 
 */
char* str_delim_max (const char *string, const char *delimeter)
{
    if (string == NULL || delimeter == NULL)
    {
        fprintf(stderr, "Пустой указатель\n");
        return NULL;
    }
    char *buffer = calloc(strlen(string) + 1, sizeof(char));
    if (buffer == NULL)
    {
        fprintf(stderr, "Ошибка выделения дин. памяти\n");
        return NULL;
    }
    const char *pt_string = string;
    size_t maxlen = 0;
    while (*pt_string)
    {
        pt_string += strspn(pt_string, delimeter);
        if (!(*pt_string))
        {
            break;
        }
        size_t len = strcspn(pt_string, delimeter);
        if (len > maxlen)
        {
            maxlen = len;
            memcpy(buffer, pt_string, len);
            buffer[len] = '\0';
        }
        pt_string += len;
    }
    return buffer;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    if (read_file_g10(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    char *delim = str_delim_max(buffer, " ");
    if (delim == NULL)
    {
        return EXIT_FAILURE;
    }
    write_to_file_g10(FILE_OUTPUT, delim);
    free(delim);
    return EXIT_SUCCESS;
}
#endif