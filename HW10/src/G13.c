/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G13: В файле input.txt записан полный адрес файла (возможно, без расширения). 
 * Необходимо изменить его расширение на ".html" и записать результат в файл output.txt. 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G13/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

char* read_file_g13 (const char *file_name, char buffer[], const int size)
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

int write_to_file_g13 (const char *file_name, const char *string)
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

const char* change_extension (char src_path[], size_t size_path, const char *extension)
{
    #define SLASH '/'
  
    if (src_path == NULL || extension == NULL)
    {
        fprintf(stderr, "Пустой указатель path\n");
        return NULL; 
    }
    char *last_slash = strrchr(src_path, SLASH);
    if (last_slash == NULL)
    {
        return NULL;
    }

    /* Сколько осталось места */
    size_t used = strlen(src_path);
    size_t remaining = size_path - used;

    if (*(last_slash + 1))
    {
        char *filename = last_slash + 1;
        char *extension_pos = strrchr(filename, '.');
        if (extension_pos == NULL)
        {
            strcat(filename, extension);
        }
        else
        {
            size_t old_ext_len = strlen(filename);
            size_t new_ext_len = strlen(extension);
            
            if (new_ext_len > old_ext_len + remaining)
            {
                fprintf(stderr, "change_extension() переполнение буфера\n");
                return NULL;
            }
            memset(extension_pos, 0, old_ext_len);
            strncat(src_path, extension, remaining);
        }
    }
    return src_path;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char path[MAX_LEN] = {'\0'};
    if (read_file_g13(FILE_INPUT, path, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    const char *extpath = change_extension(path, MAX_LEN, ".html");
    if (extpath == NULL)
    {
        exit(EXIT_FAILURE);
    }
    write_to_file_g13(FILE_OUTPUT, extpath);
    return EXIT_SUCCESS;
}
#endif