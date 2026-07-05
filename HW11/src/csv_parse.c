#include "csv_parse.h"
#include "temp_api.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

/**
 * private функция для записи сообщения об ошибке в структуру CsvParseStatus
 * @param status Указатель на структуру CsvParseStatus, в которую будет записано сообщение об ошибке
 * @param error_msg Сообщение об ошибке
 * @param error_code Код ошибки из перечисления CsvParseErrors
 * @param row_index Индекс строки, в которой произошла ошибка
 * @return Возвращает указатель на структуру CsvParseStatus, в которую было записано сообщение об ошибке
 */
static struct CsvParseStatus *write_error_msg 
(struct CsvParseStatus *status, const char *error_msg, CsvParseErrors error_code, size_t row_index)
{
    if (status)
    {
        size_t msg_length = strlen(error_msg) < LEN_ERR_MSG ? strlen(error_msg) : LEN_ERR_MSG - 1;

        strncpy(status->error_message, error_msg, msg_length);
        status->error_message[msg_length] = '\0';
        status->error_code = error_code;
        status->error_row = row_index;
    }
    return status;
}

/**
 * private функция для обработки разделителя fsm
 */
static void except_delimeter (struct CsvParseContext *context)
{
    context->buffer[context->length] = '\0';
    context->length = 0;
    context->current_column++;
}

/**
 * private функция для обработки перевода строки fsm
 * @attention: После вызова этой функции необходимо вызвать write_to_stats_array() для записи данных
 *             в массив TemperatureStats, иначе данные текущей строки будут потеряны.
 */
static void except_newline (struct CsvParseContext *context)
{
    context->buffer[context->length] = '\0';
    context->length = 0;
    context->current_row++;
    context->current_column = 0;
}

/**
 * Прогрессбар во время парсинга файла
 */
static void print_progress_bar(long current, long total)
{
    if (total <= 0)
    {
        return;
    }

    // процент выполнения
    int percentage = (int)((current * 100) / total);
    
    // ширина полосы индикатора (количество символов '=' внутри [ ])
    const int bar_width = 30;
    int progress_position = (int)((current * bar_width) / total);

    printf(BLUE"\rПарсинг: ["RESET);
    for (int i = 0; i < bar_width; ++i)
    {
        if (i < progress_position)
        {
            printf(BLUE"="RESET);
        }
        else if (i == progress_position)
        {
            printf(BLUE">"RESET);
        }
        else
            printf(" ");
    }
    printf(BLUE"] %3d%%"RESET, percentage);
    fflush(stdout); 
}

/**
 * private функция для записи данных в массив TemperatureStats
 * с проверкой на корректность данных.
 * @attention: Если данные некорректны, в структуру CsvParseStatus будет записано сообщение об ошибке.
 * @return Возвращает 1, если данные корректны и успешно записаны, иначе возвращает 0.
 */
static int write_to_stats_array 
(struct CsvParseContext *context, struct TemperatureStats *stats_array, size_t max_size)
{
    if (context->current_row >= max_size)
    {
        fprintf(stderr, RED "Ошибка: Индекс строки %zu превышает размер массива %zu\n" RESET, 
                context->current_row, max_size);
        return 0;
    }
    context->buffer[context->length] = '\0';
    struct TemperatureStats *current = &stats_array[context->current_row];
    struct CsvParseStatus *current_status = &context->status[context->current_row];
    char *endptr;
    int is_valid = 1;

    if (context->length == 0)
    {
        is_valid = 0;
    }
    else if (context->current_column >= 0 && context->current_column <= 5)
    {
        long value = strtol(context->buffer, &endptr, 10);
        
        if (*endptr != '\0' && !isspace((unsigned char)*endptr))
        {
            is_valid = 0;
        }
        
        if (context->current_column == 0)
        {
            current->year = is_valid ? (uint16_t)value : (uint16_t)0;
        }
        else if (context->current_column == 1)
        {
            current->month = is_valid ? (uint16_t)value : (uint16_t)13;
        }
        else if (context->current_column == 2)
        {
            current->day = is_valid ? (uint16_t)value : (uint16_t)33;
        }
        else if (context->current_column == 3)
        {
            current->hours = is_valid ? (uint16_t)value : (uint16_t)25;
        }
        else if (context->current_column == 4)
        {
            current->minutes = is_valid ? (uint16_t)value : (uint16_t)61;
        }
        else if (context->current_column == 5)
        {
            current->temperature = is_valid ? (int16_t)value : -999;
        }
    }
    else
    {
        is_valid = 0;
    }

    if (!is_valid)
    {
        char error_format_msg[LEN_ERR_MSG];
        snprintf(
            error_format_msg, 
            LEN_ERR_MSG, 
            RED "Ошибка валидации: Строка %zu, Колонка %d. Неверный формат данных: \"%s\"\n" RESET, 
            context->current_row + 1, 
            context->current_column + 1,
            context->buffer
        );
        write_error_msg(current_status, error_format_msg, CSV_PARSE_INVALID_FORMAT, context->current_row);
    }
    
    return is_valid;
}

/**
 * private парсер CSV-файла
 * @attention: Функция не выделяет память под массив CsvParseStatus, он должен быть выделен до вызова функции.
 * @return Возвращает указатель на массив структур CsvParseStatus, содержащий информацию об ошибках.
 */
static struct CsvParseStatus *_parse 
(FILE *file, struct TemperatureStats *stats_array, size_t size, struct CsvParseContext *context, long file_size)
{
    struct CsvParseStatus *pstatus = context->status;
    if (!file || !stats_array || !pstatus)
    {
        return pstatus;
    }

    int character;
        
    while ((character = fgetc(file)) != EOF)
    {
        if (context->current_row >= size)
        {
            write_error_msg(&pstatus[size - 1], 
                            RED "Ошибка: Количество строк в файле превышает размер массива\n" RESET, 
                            CSV_PARSE_OUT_OF_BOUNDS, 
                            context->current_row);
            return pstatus;
        }

        if (context->length >= MAX_FIELD_SIZE - 1 &&
            character != *context->delimiter && 
            character != '\n' &&
            character != '\r')
        {
            write_error_msg(&pstatus[context->current_row], 
                            RED "Ошибка: Размер поля превышает MAX_FIELD_SIZE\n" RESET, 
                            CSV_PARSE_OUT_OF_BOUNDS, 
                            context->current_row);
            return pstatus;
        }
        if (character == '\r')
        {
            continue;
        }
        if (character == *context->delimiter)
        {
            write_to_stats_array(context, stats_array, size);
            except_delimeter(context);
        }
        else if (character == '\n')
        {
            write_to_stats_array(context, stats_array, size);
            except_newline(context);
            if (context->current_row % 10000 == 0) 
            {
                print_progress_bar(ftell(file), file_size);
            }
        }
        else
        {
            context->buffer[context->length++] = (char)character;
        }
    }

    // Обработка последней строки без завершающего перевода строки \n
    if (context->length > 0 && context->current_row < size)
    {
        context->buffer[context->length] = '\0';
        write_to_stats_array(context, stats_array, size);
        context->current_row++;
    }

    print_progress_bar(file_size, file_size);
    printf("\n");

    return pstatus;
}

/* ---------------------------------- public api ------------------------------------------- */

/**
 * public функция для подсчета количества ошибок в массиве CsvParseStatus
 * @param status_array Указатель на массив структур CsvParseStatus
 * @param size Размер массива
 * @return Возвращает количество ошибок в массиве CsvParseStatus
 */
size_t errors_numbers(struct CsvParseStatus *status_array, size_t size)
{
    size_t error_count = 0;
    for (size_t i = 0; i < size; ++i)
    {
        if (strlen(status_array[i].error_message) > 0)
        {
            error_count++;
        }
    }
    return error_count;
}

/**
 * public функция для вывода сообщений об ошибках из массива CsvParseStatus
 * @param status_array Указатель на массив структур CsvParseStatus
 * @param size Размер массива
 * @return void
 */
void print_errors(struct CsvParseStatus *status_array, size_t size)
{
    for (size_t i = 0; i < size; ++i)
    {
        if (strlen(status_array[i].error_message) > 0)
        {
            fprintf(stderr, MAGENTA "Строка %zu: %s" RESET, status_array[i].error_row + 1, status_array[i].error_message);
        }
    }
}

/**
 * public парсер CSV-файла
 * @attention ВНИМАНИЕ!: Функция выделяет память под массив размером size.
 * @return Возвращает указатель на массив структур CsvParseStatus, содержащий информацию об ошибках.
 *         Его необходимо освободить после использования с помощью free().
 * @param filename Указатель на имя CSV-файла
 * @param stats_array Указатель на массив структур TemperatureStats, куда будут записаны данные
 * @param size Размер массива stats_array
 */
struct CsvParseStatus *parse_csv 
(const char *filename, struct TemperatureStats *stats_array, size_t size)
{
    struct CsvParseStatus *status_array = calloc(size, sizeof(struct CsvParseStatus));
    if (status_array == NULL)
    {
        fprintf(stderr, RED "Ошибка выделения памяти для статус-массива парсинга\n" RESET);
        return NULL;
    }

    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        fprintf(stderr, RED "Ошибка: Не удалось открыть файл %s\n" RESET, filename);
        free(status_array);
        return NULL;
    }
    // буфер для оптимизации ввода-вывода
    setvbuf(file, NULL, _IOFBF, 65536);

    // Размер файла
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);

    struct CsvParseContext context = {
                                        .status = status_array,
                                        .delimiter = ";",
                                        .length = 0,
                                        .current_row = 0,
                                        .current_column = 0,
                                    };
                                    
    _parse(file, stats_array, size, &context, file_size);
    fclose(file);

    return status_array;
}
