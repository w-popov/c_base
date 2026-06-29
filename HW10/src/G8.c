/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G8: В файле input.txt дана строка, не более 1000 символов, содержащая буквы, 
 * целые числа и иные символы. Требуется все числа, которые встречаются в строке, 
 * поместить в отдельный целочисленный массив. 
 * Например, если дана строка "data 48 call 9 read13 blank0a", то в массиве числа 48, 9, 13 и 0. 
 * Вывести массив по возрастанию в файл output.txt. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/*
  CATALOG_FILES_DEF_HW10 - полный путь к 
  каталогу файлов .txt у меня. 
  Определяется в HW10/src/CMakeLists.txt
*/
#ifdef CATALOG_FILES_DEF_HW10
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G8/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 1001

typedef int (*Compare)(int, int);

int asc (int a, int b)
{
    return a < b;
}

int *sort_X_array (int array[], const int size_array, Compare cmp)
{
    if (array == NULL)
    {
        fprintf(stderr, "Пустой указатель array[]\n");
        return NULL; 
    }
    for (int i = 1; i < size_array; ++i)
    {
        int key = array[i];
        int j = i - 1;
        while (j >= 0 && cmp(key, array[j]))
        {
            array[j + 1] = array[j];
            j--;
        }
        array[j + 1] = key;
    }
    return array;
}

char* read_file_g8 (const char *file_name, char buffer[], const int size)
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

int write_to_file_g8 (const char *file_name, int *numbers, const int size)
{
    FILE *file = NULL;

    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return 0;
    }
    
    for (int i = 0; i < size; ++i)
    {
        if (numbers[i] != INT32_MIN)
        {
            fprintf(file, "%d ", numbers[i]);
        }
    }
    fputc('\n', file);
    fclose(file);
    return 1;
}

int *nums_array (const int size)
{
    int *numbers = malloc(size * sizeof(int));
    for (int i = 0; i < size; ++i)
    {
        numbers[i] = INT32_MIN;
    }
    if (numbers == NULL)
    {
        fprintf(stderr, "Ошибка выделения int *numbers\n");
        return NULL;
    }
    return numbers;
}

int *extract_numbers(const char buffer[], int *numbers)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL; 
    }
    const char *ptbuf = buffer;
    while (*ptbuf != '\0')
    {
        char *endptr;
        long num = strtol(ptbuf, &endptr, 10);
        if (ptbuf != endptr)
        {
            *numbers++ = (int)num;
            ptbuf = endptr;
        }
        else
        {
            ptbuf++;
        }
    }
    return numbers;
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    int *ext_numbers = nums_array(MAX_LEN);
    if (ext_numbers == NULL)
    {
        return EXIT_FAILURE;
    }
    if (read_file_g8(FILE_INPUT, buffer, MAX_LEN) == NULL)
    {
        return EXIT_FAILURE;
    }
    extract_numbers(buffer, ext_numbers);
    sort_X_array(ext_numbers, MAX_LEN, asc);
    write_to_file_g8(FILE_OUTPUT, ext_numbers, MAX_LEN);

    free(ext_numbers);
    return EXIT_SUCCESS;
}
#endif