#include "csv_parse.h"
#include "temp_api.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

static void except_delimeter (struct CsvParseContext *context)
{
    context->buffer[context->length] = '\0';
    context->length = 0;
    context->current_column++;
}

static void except_newline (struct CsvParseContext *context)
{
    context->buffer[context->length] = '\0';
    context->length = 0;
    context->current_row++;
    context->current_column = 0;
}

/**
 * private функция для записи данных в массив TemperatureStats
 * с проверкой на корректность данных
 */
static void write_to_stats_array (struct CsvParseContext *context,
                                  struct TemperatureStats *stats_array,
                                  size_t max_size)
{
    if (context->current_row >= max_size)
    {
        context->error_code = CSV_PARSE_OUT_OF_BOUNDS;
        return;
    }

    struct TemperatureStats *current = &stats_array[context->current_row];
    char *endptr;
    int is_valid = 1;

    if (context->length == 0)
    {
        is_valid = 0;
    }
    else if (context->current_column >= 0 && context->current_column <= 4)
    {
        long val = strtol(context->buffer, &endptr, 10);
        if (*endptr != '\0')
        {
            is_valid = 0;
        }
        
        if (context->current_column == 0)
        {
            current->year = is_valid ? (uint16_t)val : (uint16_t)9999;
        }
        else if (context->current_column == 1)
        {
            current->month = is_valid ? (uint16_t)val : (uint16_t)99;
        }
        else if (context->current_column == 2)
        {
            current->day = is_valid ? (uint16_t)val : (uint16_t)99;
        }
        else if (context->current_column == 3)
        {
            current->hours = is_valid ? (uint16_t)val : (uint16_t)99;
        }
        else if (context->current_column == 4)
        {
            current->minutes = is_valid ? (uint16_t)val : (uint16_t)99;
        }
    }
    else if (context->current_column == 5)
    {
        float val = strtof(context->buffer, &endptr);

        if (*endptr != '\0')
        {
            is_valid = 0;
        }
        current->temperature = is_valid ? val : -999.0f;
    }

    if (!is_valid)
    {
        printf(RED "Ошибка валидации: Строка %zu, Колонка %d. Неверный формат данных: %s\n" RESET,
               context->current_row + 1, 
               context->current_column + 1,
               context->buffer
            );

        context->error_code = CSV_PARSE_INVALID_FORMAT;
    }
    memset(context->buffer, '\0', MAX_FIELD_SIZE);
}

/**
 * private парсер CSV-файла
 */
static struct TemperatureStats *_parse (FILE *file,
                                        struct TemperatureStats *stats_array,
                                        size_t size,
                                        struct CsvParseContext *context)
{
    if (!file || !stats_array || size == 0 || !context)
    {
        if (context)
        {
            context->error_code = (!file || !stats_array)
                                      ? CSV_PARSE_NULL_POINTER_ARGUMENT
                                      : CSV_PARSE_ZERO_SIZE_ARGUMENT;
        }
        return NULL;
    }

    int character;
    States state = STATE_START_FIELD;

    while ((character = fgetc(file)) != EOF)
    {
        if (context->length >= MAX_FIELD_SIZE - 1 &&
            character != *context->delimiter && character != '\n' &&
            character != '\r')
        {
            context->error_code = CSV_PARSE_INVALID_FORMAT;
            return NULL;
        }

        if (character == '\r')
        {
            continue;
        }

        switch (state)
        {
        case STATE_START_FIELD:
        case STATE_READING_FIELD:
        case STATE_EXPECT_NEXT_FIELD:
        case STATE_EXPECT_NEXT_ROW:
            if (character == *context->delimiter)
            {
                write_to_stats_array(context, stats_array, size);
                except_delimeter(context);
                state = STATE_EXPECT_NEXT_FIELD;
            }
            else if (character == '\n')
            {
                write_to_stats_array(context, stats_array, size);
                except_newline(context);
                state = STATE_EXPECT_NEXT_ROW;
            }
            else
            {
                context->buffer[context->length++] = (char)character;
                state = STATE_READING_FIELD;
            }
            break;

        default:
            break;
        }
    }

    if (context->length > 0)
    {
        write_to_stats_array(context, stats_array, size);
        context->current_row++;
    }

    return stats_array;
}

/**
 * public парсер CSV-файла
 */
struct TemperatureStats *parse_csv (const char *filename,
                                    struct TemperatureStats *stats_array,
                                    size_t size)
{
    struct CsvParseContext context = {.delimiter = ";",
                                      .length = 0,
                                      .current_row = 0,
                                      .error_row = -1,
                                      .error_column = -1,
                                      .current_column = 0,
                                      .error_code = CSV_PARSE_SUCCESS};

    char c_buffer[65536];
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        return NULL;
    }
    setvbuf(file, c_buffer, _IOFBF, sizeof(c_buffer));

    struct TemperatureStats *result = _parse(file, stats_array, size, &context);
    fclose(file);
    return result;
}
