#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "parser_csv.h"

struct SVector* svector_init (struct SVector *vec, size_t size_item, size_t cap)
{
    vec->size = 0;
    vec->capacity = cap > 0 ? cap : INITIAL_CAPACITY_SVECTOR;
    vec->size_item = size_item;
    vec->data = malloc(vec->capacity * vec->size_item);
    if (vec->data == NULL)
    {
        vec->capacity = 0;
        return NULL;
    }
    return vec;
}

void* svector_push (struct SVector *vec, void *item)
{
    if (vec->size >= vec->capacity)
    {
        vec->capacity *= REALLOC_CAPACITY_SVECTOR;
        void *new_data = realloc(vec->data, vec->capacity * vec->size_item);
        if (new_data == NULL)
        {
            return NULL;
        }
        vec->data = new_data;
    }
    char *t = (char*)vec->data + (vec->size * vec->size_item);
    memcpy(t, item, vec->size_item);
    vec->size++;
    return vec->data;
}

void* svector_get (struct SVector *vec, size_t index)
{
    if (index >= vec->size)
    {
        return NULL;
    }
    return (char*)vec->data + (index * vec->size_item);
}

void svector_free (struct SVector *vec)
{
    free(vec->data);
    vec->data = NULL;
    vec->size = 0;
    vec->capacity = 0;
    vec->size_item = 0;
}


int get_char_from_file(void *stream)
{
    return fgetc((FILE *)stream);
}


int64_t get_pos_from_file(void *stream)
{
    return (int64_t)ftell((FILE *)stream);
}


int get_char_from_string(void *stream)
{
    char **str_ptr = (char **)stream;
    if (str_ptr == NULL || *str_ptr == NULL)
    {
        return EOF;
    }
    if (**str_ptr == '\0')
    {
        return EOF;
    } 
    
    return (unsigned char)(*(*str_ptr)++); 
}


FILE *open_file (const char *filename, long *fsize)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL || fsize == NULL)
    {
        return NULL;
    }
    
    // Размер файла
    if (fseek(file, 0, SEEK_END) != 0) 
    {
        fclose(file);
        return NULL;
    }
    long file_size = ftell(file);
    if (fseek(file, 0, SEEK_SET) != 0) 
    {
        fclose(file);
        return NULL;
    }
    
    // буфер для оптимизации ввода-вывода
    setvbuf(file, NULL, _IOFBF, 65536);
    *fsize = file_size;
    return file;
}

/**
 * private функция для обработки разделителя
 */
static void except_delimeter (struct ContextParser *context)
{
    context->buffer[context->length_field] = '\0';
    context->length_field = 0;
    context->current_column++;
}

/**
 * private функция для обработки перевода строки
 */
static void except_newline (struct ContextParser *context)
{
    context->buffer[context->length_field] = '\0';
    context->length_field = 0;
    context->current_row++;
    context->current_column = 0;
    memset(context->buffer, 0, MAX_FIELD_SIZE);
}


struct ContextParser *parse_csv (struct ContextParser *context, struct ParseSource *source)
{
    if (!context || !source)
    {
        return NULL;
    }

    int character;
    int is_field_overflow = 0; 

    while ((character = source->get_char(source->stream)) != EOF)
    {
        if (character == '\r')
        {
            continue;
        }
        if (!is_field_overflow && context->length_field >= MAX_FIELD_SIZE - 1)
        {
            if (character != *context->delimiter && character != '\n')
            {
                is_field_overflow = 1;
                
                struct ErrorParse error;
                snprintf(
                    error.error_message, 
                    LEN_ERR_MSG, 
                    "Ошибка валидации: Строка %zu, Колонка %d. Размер поля превышает MAX_FIELD_SIZE\n", 
                    context->current_row + 1, 
                    context->current_column + 1
                );
                error.error_column = context->current_column + 1;
                error.error_row = context->current_row + 1;
                svector_push(context->errors_parse, &error);
            }
        }
        if (character == *context->delimiter)
        {
            if (context->current_column >= context->nums_field)
            {
                struct ErrorParse error;
                snprintf(
                    error.error_message, 
                    LEN_ERR_MSG, 
                    "Ошибка структуры: Строка %zu содержит избыточные поля (обнаружена лишняя колонка %d)\n", 
                    context->current_row + 1, 
                    context->current_column + 1
                );
                error.error_column = context->current_column + 1;
                error.error_row = context->current_row + 1;
                svector_push(context->errors_parse, &error);
            }
            else if (!is_field_overflow)
            {
                context->buffer[context->length_field] = '\0';
                context->clb_write_to_arr(context);
            }
            
            except_delimeter(context);
            is_field_overflow = 0; // Сброс флага для следующего поля
        }
        else if (character == '\n')
        {
            if (context->current_column <= context->nums_field && !is_field_overflow)
            {
                context->buffer[context->length_field] = '\0';
                context->clb_write_to_arr(context);
            }

            // проверка на нехватку полей
            if (context->current_column < context->nums_field)
            {
                struct ErrorParse error;
                snprintf(
                    error.error_message, 
                    LEN_ERR_MSG, 
                    "Ошибка структуры: Строка %zu содержит меньше полей (%d из %d)\n", 
                    context->current_row + 1, 
                    context->current_column + 1,
                    context->nums_field + 1
                );
                error.error_column = context->current_column + 1;
                error.error_row = context->current_row + 1;
                svector_push(context->errors_parse, &error);
            }
            // проверка на переизбыток полей
            else if (context->current_column > context->nums_field)
            {
                struct ErrorParse error;
                snprintf(
                    error.error_message, 
                    LEN_ERR_MSG, 
                    "Ошибка структуры: Строка %zu содержит больше полей (%d из %d)\n", 
                    context->current_row + 1, 
                    context->current_column + 1,
                    context->nums_field + 1
                );
                error.error_column = context->current_column + 1;
                error.error_row = context->current_row + 1;
                svector_push(context->errors_parse, &error);
            }
            
            except_newline(context);
            is_field_overflow = 0; 
            
            if (context->current_row % 10000 == 0 && source->get_pos && context->clb_progress) 
            {
                context->clb_progress(source->get_pos(source->stream), context->file_size);
            }
        }
        else
        {
            // Накапливать символы только для разрешенных колонок и без переполнения буфера
            if (!is_field_overflow && context->current_column <= context->nums_field)
            {
                context->buffer[context->length_field++] = (char)character;
            }
        }
    }
    
    // финал файла (если файл завершился без переноса строки)
    if (context->length_field > 0 && !is_field_overflow)
    {
        if (context->current_column <= context->nums_field)
        {
            context->buffer[context->length_field] = '\0';
            context->clb_write_to_arr(context);
        }
        
        // проверка на нехватку или избыток полей
        if (context->current_column < context->nums_field)
        {
            struct ErrorParse error;
            snprintf(
                error.error_message, 
                LEN_ERR_MSG, 
                "Ошибка структуры: Финальная строка %zu содержит меньше полей (%d из %d)\n", 
                context->current_row + 1, 
                context->current_column + 1,
                context->nums_field + 1
            );
            error.error_column = context->current_column + 1;
            error.error_row = context->current_row + 1;
            svector_push(context->errors_parse, &error);
        }
        else if (context->current_column > context->nums_field)
        {
            struct ErrorParse error;
            snprintf(
                error.error_message, 
                LEN_ERR_MSG, 
                "Ошибка: Финальная строка %zu содержит больше полей (%d из %d)\n", 
                context->current_row + 1, 
                context->current_column + 1,
                context->nums_field + 1
            );
            error.error_column = context->current_column + 1;
            error.error_row = context->current_row + 1;
            svector_push(context->errors_parse, &error);
        }
        
        context->current_row++;
    }
    // Если файл завершился пустым полем после ';'
    else if (context->current_column > 0)
    {
        if (context->current_column < context->nums_field)
        {
            struct ErrorParse error;
            snprintf(
                error.error_message, 
                LEN_ERR_MSG, 
                "Ошибка: Строка %zu оборвалась, обнаружено меньше полей (%d из %d)\n", 
                context->current_row + 1, 
                context->current_column + 1,
                context->nums_field + 1
            );
            error.error_column = context->current_column + 1;
            error.error_row = context->current_row + 1;
            svector_push(context->errors_parse, &error);
        }
        else if (context->current_column > context->nums_field)
        {
            struct ErrorParse error;
            snprintf(
                error.error_message, 
                LEN_ERR_MSG, 
                "Ошибка: Строка %zu содержит больше полей (%d из %d)\n", 
                context->current_row + 1, 
                context->current_column + 1,
                context->nums_field + 1
            );
            error.error_column = context->current_column + 1;
            error.error_row = context->current_row + 1;
            svector_push(context->errors_parse, &error);
        }
    }
    
    if (source->get_pos && context->clb_progress)
    {
        context->clb_progress(source->get_pos(source->stream), context->file_size);
    }
    printf("\n");

    return context;
}

void show_errors(struct SVector *errors, size_t rows)
{
    size_t errnums = errors->size;
    printf("Валидных строк:...%zu\n", rows);
    if (errnums > 0)
    {
        printf("Ошибок:...........%s%zu%s\n",RED, errnums, RESET);
        for (size_t i = 0; i < errnums; ++i)
        {
            struct ErrorParse *item_err = (struct ErrorParse*)svector_get(errors, i);
            printf(RED"%s"RESET, item_err->error_message);
        }
    }
    else
    {
        printf("Ошибок:...........%s%zu%s\n",GREEN, errnums, RESET);
    }
}















