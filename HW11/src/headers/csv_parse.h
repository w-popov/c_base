#ifndef _CSV_PARSE_H_
#define _CSV_PARSE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#define MAX_FIELD_SIZE 64

typedef enum {
    STATE_START_FIELD,
    STATE_READING_FIELD,
    STATE_EXPECT_NEXT_FIELD,
    STATE_EXPECT_NEXT_ROW,
    STATE_END_OF_FILE,
    STATE_ERROR
} States;

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

struct CsvParseContext
{
    const char* delimiter;          // Разделитель
    size_t length;                  // Размер ячейки
    size_t current_row;             // Индекс строки
    char buffer[MAX_FIELD_SIZE];    // Буфер для хранения текущей ячейки
    int32_t error_row;              // Номер строки ошибки
    int32_t error_column;           // Номер колонки ошибки
    int32_t current_column;         // Текущая колонка
    CsvParseErrors error_code;      // Код ошибки
};

struct TemperatureStats* parse_csv (const char* filename, 
                                    struct TemperatureStats *stats_array, 
                                    size_t size
                                    );

#ifdef __cplusplus
}
#endif

#endif