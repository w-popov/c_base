/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G4: В файле input.txt даны два слова не более 100 символов каждое, 
 *     разделенные одним пробелом. Найдите только те символы слов, 
 *     которые встречаются в обоих словах только один раз. 
 *     Напечатайте их через пробел в файл output.txt в лексикографическом порядке. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/*
  CATALOG_FILES_DEF_HW10 - полный путь к 
  каталогу файлов .txt у меня. 
  Определяется в HW10/src/CMakeLists.txt
*/
#ifdef CATALOG_FILES_DEF_HW10
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G4/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH     "input.txt"
#define FILE_OUTPUT CAT_PATH    "output.txt"
#define MAX_LEN                 101
#define MAX_WORDS               2
#define ALPHBET_RANGE           123

char (*read_words (char words[][MAX_LEN], const char* file_name))[MAX_LEN]
{
    if (words == NULL)
    {
        fprintf(stderr, "Пустой указатель words[][]\n");
        return NULL; 
    }

    FILE *file = NULL;
    if ((file = fopen(file_name, "r")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", file_name);
        return NULL;
    }
    for (int i = 0; i < MAX_WORDS; ++i)
    {
        if (fscanf(file, "%s", words[i]) < 0)
        {
            break;
        }
        words[i][strlen(words[i])] = '\0';
    }
    fclose(file);
    return words;
}

/**
 * Кол-во вхождений символа в строке
 */
int char_count(const char *str, char ch)
{
    int count = 0;
    while ((str = strchr(str, ch)) != NULL)
    {
        count++;
        str++;
    }
    return count;
}

/**
 * Выделяет память!
 * Решение задания.
 */
int *letter_matches (char words[][MAX_LEN])
{
    int *matches = calloc(ALPHBET_RANGE, sizeof(int));
    if (matches == NULL)
    {
        fprintf(stderr, "Ошибка выделения int *matches\n");
        return NULL;
    }

    const char *alphabet = "abcdefghijklmnopqrstuvwxyz";
    for (size_t i = 0; i < strlen(alphabet); ++i)
    {
        if ((char_count(words[0], alphabet[i]) == 1) && (char_count(words[1], alphabet[i]) == 1))
        {
            size_t index = (unsigned)alphabet[i];
            matches[index] += 1;
        }
    }
    return matches;
}

void write_to_file_g4 (const int *matches, const char *filename)
{
    if (matches == NULL)
    {
        fprintf(stderr, "Пустой указатель matches\n");
        return; 
    }

    FILE *file = NULL;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "Ошибка открытия: %s\n", filename);
        return;
    }
    
    for (int i = 0; i < ALPHBET_RANGE; ++i)
    {
        if (matches[i] == 1)
        {
            fprintf(file, "%c ", (char)i);
        }
    }
    fprintf(file, "\n");
    
    fclose(file);
}

#ifndef TEST_DEF_HW10
int main (void)
{
    char words[MAX_WORDS][MAX_LEN] = {'\0'}; 
    if (read_words(words, FILE_INPUT) == NULL)
    {
        return EXIT_FAILURE;
    }

    int* const match = letter_matches(words);
    if (match == NULL)
    {
        return EXIT_FAILURE;
    }

    write_to_file_g4(match, FILE_OUTPUT);
    free(match);
    return EXIT_SUCCESS;
}
#endif