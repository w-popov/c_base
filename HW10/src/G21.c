/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G21: Во входном файле input.txt записано некоторое количество 
 * символов * (камней). Необходимо построить равносторонний треугольник 
 * используя все символы * и символ пробел, записать ответ в выходной файл output.txt. 
 * Между соседними символами * строго один пробел. Если треугольник невозможно составить, 
 * используя все камни, то необходимо записать единственное слово NO в файл output.txt. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
  CATALOG_FILES_DEF_HW10 - полный путь к
  каталогу файлов .txt у меня.
  Определяется в HW10/src/CMakeLists.txt
*/
#ifdef CATALOG_FILES_DEF_HW10
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G21/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 4001

char *read_file_g21 (const char *file_name, char buffer[], const int size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL;
    }

    // "%4000[^\1] - читать всё не более 4000 символов"
    char format[18];
    snprintf(format, sizeof(format), "%%%d[^\1]", size - 1);

    FILE *file = NULL;

    if ((file = fopen(file_name, "r")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return NULL;
    }
    fscanf(file, format, buffer);
    fclose(file);

    return buffer;
}

unsigned is_triangle (int star_count)
{
    if (!star_count)
    {
        return 0;
    }
    unsigned D = 1 + 8 * (unsigned)star_count;
    unsigned sqrt_D = (unsigned)round(sqrt((double)D));

    if (sqrt_D * sqrt_D != D || (sqrt_D - 1) % 2 != 0) 
    {
        return 0;
    }
    return sqrt_D;
}

/**
 * Чтобы составить равносторонний треугольник из всех камней (*)
 * их общее количество должно быть треугольным числом, 
 * то есть выражаться формулой N = k * (k + 1) / 2 
 * где k - количество камней в основании,
 *     N - количество звезд в строке.
 */
void equilateral_triangles (char string[], const char *file_name)
{
    int star_count = 0;
    const char *pstr = string;
    while (*pstr) 
    {
        if (*pstr == '*') 
        {
            star_count++;
        }
        ++pstr;
    }
    FILE *file = NULL;

    if ((file = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return;
    }
    
    unsigned sqrt_d = is_triangle(star_count);
    if (!sqrt_d)
    {
        fprintf(file, "%s\n", "NO");
        goto end;
    }

    // k — высота треугольника (количество строк)
    size_t k = (sqrt_d - 1) / 2;

    for (size_t i = 1; i <= k; ++i) 
    {
        for (size_t space = 0; space < k - i; ++space) 
        {
            fputc(' ', file);
        }
        for (size_t star = 0; star < i; ++star) 
        {
            fputc('*', file);
            if (star < i - 1) 
            {
                fputc(' ', file);
            }
        }
        fputc('\n', file);
    }

    end:
    fclose(file);
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char string[MAX_LEN] = {'\0'};
    read_file_g21(FILE_INPUT, string, MAX_LEN);
    equilateral_triangles(string, FILE_OUTPUT);
    return EXIT_SUCCESS;
}
#endif