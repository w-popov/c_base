/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G14: В файле input.txt в одной строке фамилию, имя и отчество. 
 * Сформировать файл приветствие output.txt, где останутся имя и фамилия 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G14/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 101

char* read_file_g14 (const char *file_name, char buffer[], const int size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL; 
    }

    // "%100[^\n]"
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

int write_to_file_g14 (const char *file_name, const char *string)
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
 * @attention Выделяет и освобождает память!
 */
void hello_name_patronymic (char human[], const int size, const char *filename)
{
    if (human == NULL)
    {
        fprintf(stderr, "Пустой указатель\n");
        return;
    }
    #define HUMAN_LEN size
    char *hello_name = calloc(size + 16, sizeof(char));
    char *phuman = human;
    char name[HUMAN_LEN];          // имя
    char lastname[HUMAN_LEN];      // фамилия
    char patronymic[HUMAN_LEN];    // отчество

    strcat(hello_name, "Hello, ");
    for (int i = 0; *phuman; ++i)
    {
        phuman += strspn(phuman, " ");
        if (!(*phuman))
        {
            break;
        }
        size_t len = strcspn(phuman, " ");
        switch (i)
        {
        case 0:
            memcpy(lastname, phuman, len);
            lastname[len] = '\0';
            break;
        case 1:
            memcpy(name, phuman, len);
            name[len] = '\0';
            break;
        case 2:
            memcpy(patronymic, phuman, len);
            patronymic[len] = '\0';
            break;

        default:
            break;
        }
        phuman += len;
    }
    strcat(hello_name, name);
    strcat(hello_name, " ");
    strcat(hello_name, lastname);
    strcat(hello_name, "!");

    write_to_file_g14(filename, hello_name);
    free(hello_name);
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char human[MAX_LEN] = {'\0'};
    if (read_file_g14(FILE_INPUT, human, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    hello_name_patronymic(human, MAX_LEN, FILE_OUTPUT);
    return EXIT_SUCCESS;
}
#endif