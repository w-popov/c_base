/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G2: Считать число N из файла input.txt. Сформировать строку из N символов.
 *     N четное число, не превосходящее 26. На четных позициях должны находится 
 *     четные цифры в порядке возрастания, кроме 0, на нечетных позициях - 
 *     заглавные буквы в порядке следования в английском алфавите. 
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G2/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 64
#ifndef TEST_DEF_HW10

int read_file_g2 (const char *filename)
{
    FILE *file = NULL;
    if ( (file = fopen(filename, "r")) == NULL )
    {
        fprintf(stderr, "Ошибка открытия: %s\n", filename);
        return 0;
    }
    int N = 0;
    if (fscanf(file, "%d", &N) < 0)
    {
        fprintf(stderr, "Ошибка чтения или пустой файл: %s\n", filename);
        fclose(file);
        return 0;
    }
    fclose(file);
    return N;
}

const char* line_and_numbers (int number, char buffer[], const int size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL;
    }
    if (number > 26 || number == 0)
    {
        fprintf(stderr, "Число вне диапазона 1..26\n");
        return NULL;
    }
   
    int i = 1, count_even = 0;
    for (int litera = 65; i <= number && litera < 91 && i < size; ++i)
    {
        if (!(i % 2))
        {
            ++count_even;
            buffer[i-1] = (2 * ((count_even - 1) % 4) + 2) + '0';
        }
        else 
        {
            buffer[i-1] = (char)(litera++);
        }
    }
    buffer[i] = '\0';
    return buffer;
}

void write_to_file_g2 (const char *filename, const char buffer[])
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return; 
    }

    FILE *file = NULL;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", filename);
        return;
    }
    fprintf(file, "%s\n", buffer);
    
    fclose(file);
}

int main (void)
{
    char buffer[MAX_LEN] = {'\0'};
    int number = read_file_g2(FILE_INPUT);
    if(line_and_numbers(number, buffer, MAX_LEN) == NULL)
    {
        fprintf(stderr, "line_and_number() вернула NULL\n");
        return EXIT_FAILURE; 
    }
    write_to_file_g2(FILE_OUTPUT, buffer);
    return EXIT_SUCCESS;
}
#endif