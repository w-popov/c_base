/**
 * ДЗ-10. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * G22: Известный алгоритм Soundex...
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
#define CAT_PATH CATALOG_FILES_DEF_HW10 "G22/"
#else
#define CAT_PATH ""
#endif

#define FILE_INPUT CAT_PATH "input.txt"
#define FILE_OUTPUT CAT_PATH "output.txt"
#define MAX_LEN 21

typedef int (*State)(char *, char *, int *);

enum {
    STOP,
    BEGIN,
    FIRST_CH,
    VOWEL_SOUND,
    CONSON_1,
    CONSON_2,
    CONSON_3,
    CONSON_4,
    CONSON_5,
    CONSON_6,
    END
};

int what (char c)
{
    c = tolower(c);

    enum { ROW = 7, COL = 8 };
    char sounds[ROW][COL] = {
        {'a', 'e', 'h', 'i', 'o', 'u', 'w', 'y'},        // vowel
        {'b', 'f', 'p', 'v', '\0', '\0', '\0', '\0'},    // con1
        {'c', 'g', 'j', 'k', 'q', 's', 'x', 'z'},        // con2
        {'d', 't', '\0', '\0', '\0', '\0', '\0', '\0'},  // con3
        {'l', '\0', '\0', '\0', '\0', '\0', '\0', '\0'}, // con4
        {'m', 'n', '\0', '\0', '\0', '\0', '\0', '\0'},  // con5
        {'r', '\0', '\0', '\0', '\0', '\0', '\0', '\0'}  // con6
    };
    int row = 0;
    int is_finded = 0;
    for (; row < ROW; ++row)
    {
        is_finded = 0;
        for (int col = 0; col < COL; ++col)
        {
            if (c == sounds[row][col])
            {
                is_finded = 1;
                break;
            }
        }
        if (is_finded)
        {
            break;
        }
    }
    switch (row)
    {
    case 0:
        return VOWEL_SOUND;
    case 1:
        return CONSON_1;
    case 2:
        return CONSON_2;
    case 3:
        return CONSON_3;
    case 4:
        return CONSON_4;
    case 5:
        return CONSON_5;
    case 6:
        return CONSON_6;
    default:
        break;
    }
    return STOP;
}

int begin (char *s, char *code, int *code_len)
{
    if (s == NULL || *s == '\0')
    {
        return STOP;
    }
    *code = tolower(*s);
    *code_len = 1;
    s++;
    if (*s == '\0')
    {
        while (*code_len < 4)
        {
            code[*code_len] = '0';
            (*code_len)++;
        }
        code[4] = '\0';
        return STOP;
    }
    return what(*s);
}

int first_ch (char *s, char *code, int *code_len)
{
    (void)code;
    (void)code_len;
    if (*s == '\0')
    {
        return END;
    }
    return what(*s);
}

int vowel_sound (char *s, char *code, int *code_len)
{
    (void)code;
    (void)code_len;
    s++;
    if (*s == '\0')
    {
        return END;
    }
    return what(*s);
}

int add_digit (char *s, char *code, int *code_len, char digit)
{
    if (*s == '\0')
    {
        return END;
    }
    if (*code_len > 0 && code[*code_len - 1] == digit)
    {
        s++;

        if (*s == '\0')
        {
            return END;
        }
        return what(*s);
    }    
    if (*code_len < 4)
    {
        code[*code_len] = digit;
        (*code_len)++;
    }
    s++;

    if (*s == '\0')
    {
        return END;
    }
    
    return what(*s);
}

int conson_1 (char *s, char *code, int *code_len)
{
    if (*s == '\0')
    {
        return END;
    }
    return add_digit(s, code, code_len, '1');
}

int conson_2 (char *s, char *code, int *code_len)
{
    if (*s == '\0')
    {
        return END;
    }
    return add_digit(s, code, code_len, '2');
}

int conson_3 (char *s, char *code, int *code_len)
{
    if (*s == '\0')
    {
        return END;
    }
    return add_digit(s, code, code_len, '3');
}

int conson_4 (char *s, char *code, int *code_len)
{
    if (*s == '\0')
    {
        return END;
    }
    return add_digit(s, code, code_len, '4');
}

int conson_5 (char *s, char *code, int *code_len)
{
    if (*s == '\0')
    {
        return END;
    }
    return add_digit(s, code, code_len, '5');
}

int conson_6 (char *s, char *code, int *code_len)
{
    if (*s == '\0')
    {
        return END;
    }
    return add_digit(s, code, code_len, '6');
}

int end_func (char *s, char *code, int *code_len)
{
    (void)s;
    while (*code_len < 4)
    {
        code[*code_len] = '0';
        (*code_len)++;
    }
    code[4] = '\0';
    return STOP;
}

void soundex (char str[], char code[])
{
    State func_state = begin;
    char *string = str;
    int code_len = 0;
    int state = func_state(string, code, &code_len);

    while (state)
    {
        switch (state)
        {
        case FIRST_CH:
            func_state = first_ch;
            break;
        case VOWEL_SOUND:
            func_state = vowel_sound;
            break;
        case CONSON_1:
            func_state = conson_1;
            break;
        case CONSON_2:
            func_state = conson_2;
            break;
        case CONSON_3:
            func_state = conson_3;
            break;
        case CONSON_4:
            func_state = conson_4;
            break;
        case CONSON_5:
            func_state = conson_5;
            break;
        case CONSON_6:
            func_state = conson_6;
            break;
        case END:
            func_state = end_func;
            break;
        default:
            func_state = end_func;
            break;
        }

        state = func_state(string, code, &code_len);
        string++;
    }
}

char *read_file_g22 (const char *file_name, char buffer[], const int size)
{
    if (buffer == NULL)
    {
        fprintf(stderr, "Пустой указатель buffer[]\n");
        return NULL;
    }

    // "%20[^\n]"
    char format[18];
    snprintf(format, sizeof(format), "%%%d[^\n]", size - 1);

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

int write_to_file_g22 (const char *file_name, const char *string)
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

#ifndef TEST_DEF_HW10
int main ()
{
    char string[MAX_LEN] = {'\0'};
    char code[MAX_LEN] = {'\0'};
    read_file_g22(FILE_INPUT, string, MAX_LEN);
    soundex(string, code);
    write_to_file_g22(FILE_OUTPUT, code);

    return 0;
}
#endif