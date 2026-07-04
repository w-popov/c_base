#ifndef _CSV_PARSE_H_
#define _CSV_PARSE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdio.h>

#define MAX_FIELD_SIZE 64   // Максимальный размер поля в CSV-файле
#define LEN_ERR_MSG 256     // Максимальная длина сообщения об ошибке

typedef enum {
    CSV_PARSE_SUCCESS = 0,
    CSV_PARSE_FILE_NOT_FOUND,
    CSV_PARSE_INVALID_FORMAT,
    CSV_PARSE_MEMORY_ALLOCATION_FAILED,
    CSV_PARSE_NULL_POINTER_ARGUMENT,
    CSV_PARSE_ZERO_SIZE_ARGUMENT,
    CSV_PARSE_OUT_OF_BOUNDS,
    SCV_PARSE_UNKNOWN_ERROR
} CsvParseErrors;

struct TemperatureStats; 

struct CsvParseStatus
{
    char error_message[LEN_ERR_MSG];    // Сообщение об ошибке
    CsvParseErrors error_code;          // Код ошибки
    size_t error_row;                   // Номер строки ошибки
    size_t rows_parsed;                 // Количество успешно обработанных строк
    int16_t error_column;               // Номер колонки ошибки
};

struct CsvParseContext
{
    struct CsvParseStatus *status;      // Статус-массива парсинга
    const char* delimiter;              // Разделитель
    size_t length;                      // Размер ячейки
    size_t current_row;                 // Индекс строки
    char buffer[MAX_FIELD_SIZE];        // Буфер для хранения текущей ячейки
    int16_t current_column;             // Текущая колонка
};

struct CsvParseStatus* parse_csv 
(const char* filename, struct TemperatureStats *stats_array, size_t size);

size_t errors_numbers(struct CsvParseStatus *status_array, size_t size);

void print_errors(struct CsvParseStatus *status_array, size_t size);

#ifdef __cplusplus
}
#endif

#endif