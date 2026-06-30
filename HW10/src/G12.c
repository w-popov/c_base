/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G12: В файле input.txt дано предложение требуется разобрать 
 * его на отдельные слова. Напечатать каждое слово на отдельной строке в файл output.txt. 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G12/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

char* read_file_g12 (const char *file_name, char buffer[], const int size)
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

/**
 * @brief Разбить строку на части по разделителю и записывает в файл
 * @details strspn - Возвращает длину начального сегмента str, 
 *                   состоящего только из символов, перечисленных в accept.
 *          strcspn - Возвращает длину начального сегмента str, 
 *                    не содержащего ни одного символа из reject
 * @param string входная строка.
 * @param delimeter разделитель.
 * @param file_name имя выходного файла
 * @return void
 * 
 */
void str_delim (const char *string, const char *delimeter, const char *file_name)
{
    if (string == NULL || delimeter == NULL || file_name == NULL)
    {
        fprintf(stderr, "Пустой указатель\n");
        return;
    }
    
    const char *pt_string = string;
    FILE *file = NULL;
    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return;
    }    
    while (*pt_string)
    {
        pt_string += strspn(pt_string, delimeter);
        if (!(*pt_string))
        {
            break;
        }
        size_t len = strcspn(pt_string, delimeter);
        fprintf(file, "%.*s\n", (int)len, pt_string);
        pt_string += len;
    }
    fclose(file);
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    if (read_file_g12(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    str_delim(buffer, " ", FILE_OUTPUT);
    return EXIT_SUCCESS;
}
#endif