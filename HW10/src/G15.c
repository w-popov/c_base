/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G15: В файле input.txt дано предложение. Необходимо заменить 
 * все имена «Cao» на «Ling» и записать результат в файл output.txt. 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G15/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

char* read_file_g15 (const char *file_name, char buffer[], const int size)
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
    if (fscanf(file, format, buffer) < 0)
    {
        fprintf(stderr, "Ошибка чтения или пустой файл: %s\n", file_name);
        fclose(file);
        return NULL;
    }
    fclose(file);
    
    return buffer;
}

int write_to_file_g15 (const char *file_name, const char *string)
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
 * @attention Выделяет память!
 */
char* cao_ling (const char *string, const char *old_word, const char *new_word)
{
    if (string == NULL)
    {
        fprintf(stderr, "Пустой указатель\n");
        return NULL;
    }

    size_t new_len = strlen(new_word);
    size_t old_len = strlen(old_word);
    int i = 0;
    size_t count = 0;
    const char *tmp = string;
    while ((tmp = strstr(tmp, old_word))) 
    {
        count++;
        tmp += old_len;
    }
    /* точный размер новой строки */
    size_t orig_len = strlen(string);
    size_t final_len = orig_len + count * (new_len - old_len);
    
    char *result = calloc(final_len + 1, sizeof(char));

    if (result == NULL) 
    {
        return NULL;
    }
    while (*string)
    {
        if (strstr(string, old_word) == string)
        {
            strcpy(&result[i], new_word);
            i += new_len;
            string += old_len;
        }
        else
        {
            result[i++] = *string++;
        }
    }
    result[i] = '\0';
    return result;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char cao_string[MAX_LEN] = {'\0'};
    if (read_file_g15(FILE_INPUT, cao_string, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    char *ling_string = cao_ling(cao_string, "Cao", "Ling");
    if (ling_string == NULL)
    {
        return EXIT_FAILURE;
    }
    write_to_file_g15(FILE_OUTPUT, ling_string);
    free(ling_string);
    return EXIT_SUCCESS;
}
#endif