#ifndef _PARSER_CSV_H_
#define _PARSER_CSV_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdio.h>

// Цвета
#define RED "\033[31m"
#define RED_BOLD "\033[1;31m"
#define GREEN "\033[32m"
#define GREEN_BOLD "\033[1;32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"
#define RESET "\033[0m"

#define INITIAL_CAPACITY_SVECTOR  1024      // Начальная ёмкость SVector
#define REALLOC_CAPACITY_SVECTOR  4         // Во сколько раз увеличится размер вектора
#define MAX_FIELD_SIZE 64                   // Максимальный размер поля в CSV-файле
#define LEN_ERR_MSG 255                     // Максимальная длина сообщения об ошибке

/**
 * Динамический массив.
 */
struct SVector
{
    void *data;             // Указатель на массив
    size_t capacity;        // Ёмкость (кол-во элементов)
    size_t size;            // Текущее кол-во элементов
    size_t size_item;       // Размер одного элемента sizeof(item)
};

/**
 * Ошибки
 */
struct ErrorParse
{
    char error_message[LEN_ERR_MSG];    // Сообщение об ошибке
    size_t error_row;                   // Номер строки ошибки
    int16_t error_column;               // Номер колонки ошибки
};


/** 
 * @brief Указатель на callback-функцию для посимвольного чтения
 */
typedef int (*GetNextCharCallback)(void *source);

/**
 * @brief Указатель на callback-функцию для прогрессбара
 */
typedef int64_t (*GetPosCallback)(void *source);

struct ParseSource 
{
    void *stream;                  // Указатель на FILE* или на строку char*
    GetNextCharCallback get_char;  // Функция чтения
    GetPosCallback get_pos;        // Функция получения текущей позиции (для прогресс-бара)
};

struct ContextParser;

/**
 * @brief Указатель на callback-функцию индикатора выполнения
 */
typedef void (*CallbackProgressBar)(long current, long total);

/**
 * @brief Указатель на callback-функцию записи в массив данных
 */
typedef int (*CallbackWriteToArray)(struct ContextParser *ctx);


struct ContextParser
{
    struct SVector *errors_parse;           // Указатель на массив структур ошибок
    struct SVector *array;                  // Указатель на массив данных 
    const char* delimiter;                  // Разделитель
    CallbackProgressBar clb_progress;       // Указатель на ф-цию прогрессбара
    CallbackWriteToArray clb_write_to_arr;  // Указатель на ф-цию записи в массив данных
    size_t file_size;                       // Размер файла
    size_t current_row;                     // Индекс строки
    char buffer[MAX_FIELD_SIZE];            // Буфер для хранения текущей ячейки   
    int16_t current_column;                 // Текущая колонка
    uint16_t length_field;                  // Размер ячейки (столбца)
    uint16_t nums_field;                    // Количество столбцов (счет с нуля)
};


/**
 * @brief Посимвольное чтение из файла
 */
int get_char_from_file(void *stream);

/**
 * @brief Посивольное чтение из строки
 */
int get_char_from_string(void *stream);

/**
 * @brief Позиция в файле для прогрессбара
 */
int64_t get_pos_from_file(void *stream);

/**
 * @brief Инициализация.
 * @param *vec указатель на вектор
 * @param size_item размер элемента вектора
 * @param cap задать ёмкость, если передано 0, то ёмкость по умолчанию
 * @return указатель на вектор, или NULL 
 */
struct SVector* svector_init (struct SVector *vec, size_t size_item, size_t cap);


/**
 * @brief Добавить один элемент.
 * @param *vec указатель на вектор
 * @param *item указатель на добавляемый элемент
 * @return *void указатель на массив data, или NULL
 */
void* svector_push (struct SVector *vec, void *item);


/**
 * @brief Получить элемент по индексу.
 * @param *vec указатель на вектор
 * @param index индекс
 * @return *void указатель на элемент в векторе
 */
void* svector_get (struct SVector *vec, size_t index);


/**
 * @brief Освободить память
 * @param *vec указатель на вектор
 * @return void
 */
void svector_free (struct SVector *vec);

/**
 * @brief Открыть файл
 * @param *filename имя файла
 * @param *fsize передача по указателю размера файла
 * @return Дескриптор открытого файла
 */
FILE* open_file (const char *filename, int64_t *fsize);

/**
 * @brief Запуск парсинга .csv файла
 * @param *ctx Указатель на контекст
 * @param *source источник данных
 * @return Указатель на контекст или NULL
 */
struct ContextParser *parse_csv (struct ContextParser *ctx, struct ParseSource *source);

/**
 * @brief Вывод ошибок
 * @param *errors вектор структур ошибок
 * @param rows кол-во строк для информативности
 */
void show_errors(struct SVector *errors, size_t rows);


#ifdef __cplusplus
}
#endif

#endif