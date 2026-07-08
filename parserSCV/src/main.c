#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "parcer_csv.h"

/*
temperature_big.csv
temperature_small.csv
*/

struct TemperatureStats
{
    uint32_t year       :11;
    uint32_t month      :4;
    uint32_t day        :6;
    uint32_t hours      :5;
    uint32_t minutes    :6;
    int16_t temperature;
};

int is_valid_number(char buffer[], int is_signed)
{
    if (buffer == NULL)
    {
        return 0;
    }

    size_t i = 0;

    while (buffer[i] != '\0' && isspace((unsigned char)buffer[i]))
    {
        i++;
    }
    if (buffer[i] == '\0')
    {
        return 0;
    }
    if (buffer[i] == '-' || buffer[i] == '+')
    {
        if (!is_signed)
        {
            return 0;
        }

        i++;
        
        if (buffer[i] == '\0' || isspace((unsigned char)buffer[i]))
        {
            return 0;
        }
    }
    if (!isdigit((unsigned char)buffer[i]))
    {
        return 0;
    }
    while (buffer[i] != '\0' && isdigit((unsigned char)buffer[i]))
    {
        i++;
    }
    while (buffer[i] != '\0' && isspace((unsigned char)buffer[i]))
    {
        i++;
    }

    return (buffer[i] == '\0');
}

int validator_fields (void *this, long value, uint16_t num_field)
{
    (void)(this);
    
    switch (num_field)
    {
    case 0:
        return (uint16_t)value > 1950 && (uint16_t)value < 2040;
    case 1:
        return (uint16_t)value > 0 && (uint16_t)value < 13;
    case 2:
        return (uint16_t)value > 0 && (uint16_t)value < 32;
    case 3:
        return (uint16_t)value < 24;
    case 4:
        return (uint16_t)value < 60;
    case 5:
        return (int16_t)value > -99 && (int16_t)value < 100;
    default:
        break;
    }
    return 0;
}

/**
 * Прогрессбар во время парсинга файла
 */
void print_progress_bar(void *this, long current, long total)
{
    (void)(this);

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

int write_to_array (void *this, struct ContextParser *context)
{
    (void)this;
    context->buffer[context->length_field] = '\0';
    
    static struct TemperatureStats current;
    static int is_row_valid = 1; // Флаг валидности всей текущей строки

    if (context->current_column == 0)
    {
        current = (struct TemperatureStats){0};
        is_row_valid = 1;
    }
    context->buffer[context->length_field] = '\0';
    int is_field_valid = 1;
    char *endptr;
    
    if (context->length_field == 0)
    {
        is_field_valid = 0;
    }
    else if (context->current_column >= 0 && context->current_column <= context->nums_field)
    {
        long value = strtol(context->buffer, &endptr, 10);
        
        if (*endptr != '\0' && !isspace((unsigned char)*endptr))
        {
            is_field_valid = 0;
        }
        
        // Валидация конкретной колонки
        if (context->current_column == 0)
        {
            if (is_field_valid && is_valid_number(context->buffer, 0) && 
                validator_fields(NULL, value, context->current_column))
            {
                current.year = value;
            }
            else
            {
                is_field_valid = 0;
            }
        }
        else if (context->current_column == 1)
        {
            if (is_field_valid && is_valid_number(context->buffer, 0) && 
                validator_fields(NULL, value, context->current_column))
            {
                current.month = value;
            }
            else
            {
                is_field_valid = 0;
            }
        }
        else if (context->current_column == 2)
        {
            if (is_field_valid && is_valid_number(context->buffer, 0) && 
                validator_fields(NULL, value, context->current_column))
            {
                current.day = value;
            }
            else
            {
                is_field_valid = 0;
            }
        }
        else if (context->current_column == 3)
        {
            if (is_field_valid && is_valid_number(context->buffer, 0) && 
                validator_fields(NULL, value, context->current_column))
            {
                current.hours = value;
            }
            else
            {
                is_field_valid = 0;
            }
        }
        else if (context->current_column == 4)
        {
            if (is_field_valid && is_valid_number(context->buffer, 0) && 
                validator_fields(NULL, value, context->current_column))
            {
                current.minutes = value;
            }
            else
            {
                is_field_valid = 0;
            }
        }
        else if (context->current_column == 5)
        {
            if (is_field_valid && is_valid_number(context->buffer, 1) && 
                validator_fields(NULL, value, context->current_column))
            {
                current.temperature = (int16_t)value;
            }
            else
            {
                is_field_valid = 0;
            }
        }
    }
    else
    {
        is_field_valid = 0;
    }    
    if (!is_field_valid)
    {
        is_row_valid = 0;
        
        struct ErrorParse error;
        snprintf(
            error.error_message, 
            LEN_ERR_MSG, 
            "Ошибка валидации: Строка %zu, Колонка %d. Неверный формат данных: \"%s\"\n", 
            context->current_row + 1, 
            context->current_column + 1,
            context->buffer
        );
        error.error_column = context->current_column + 1;
        error.error_row = context->current_row + 1;
        svector_push(context->errors_parse, &error);
    }
    if (context->current_column == context->nums_field)
    {
        if (is_row_valid)
        {
            svector_push(context->array, &current);
        }
        
        // сброс для след. строки
        current = (struct TemperatureStats){0}; 
        is_row_valid = 1;                       
    }
    
    return is_field_valid;
}



int main(void)
{
    const char *filename = "temperature_big.csv";
    long fsize = 0;
    FILE *file = open_file(filename, &fsize);

    struct SVector t_array, err_array;

    struct SVector *array = svector_init(&t_array, sizeof(struct TemperatureStats), 0);
    struct SVector *errors_array = svector_init(&err_array, sizeof(struct ErrorParse), 0);
    CallbackProgressBar progressbar = print_progress_bar;
    CallbackWriteToArray write = write_to_array;

    struct ContextParser cntx = {
        .array = array,
        .delimiter = ";",
        .file_size = fsize,
        .errors_parse = errors_array,
        .length_field = 0,
        .current_column = 0,
        .current_row = 0,
        .nums_field = 5,
        .clb_progress = progressbar,
        .clb_write_to_arr = write
    };

    
    struct ContextParser *result = parse_csv(&cntx, file);
    printf("Размер: %zu cap: %zu\n", result->array->size, result->array->capacity);

    for (size_t i = 0; i < result->errors_parse->size; ++i)
    {
        struct ErrorParse *item_err = (struct ErrorParse*)svector_get(result->errors_parse, i);
        printf("%s", item_err->error_message);
    }
    printf("---------------------%zu-------------------------------\n", result->array->size);

    
    setvbuf(stdout, NULL, _IOFBF, 65536);

    struct TemperatureStats *tarr = (struct TemperatureStats*)result->array->data;
    size_t total_items = result->array->size;

    printf("---------------------%zu-------------------------------\n", total_items);
    for (size_t i = 0; i < total_items; ++i)
    {
        printf("%6zu]  Year: %4u,  Month: %2u,  Day: %2u,  Hours: %2u,  Minutes: %2u,  Temperature: %2d\n",
            i + 1, tarr[i].year, tarr[i].month, tarr[i].day, tarr[i].hours, tarr[i].minutes, tarr[i].temperature);
    }
    fflush(stdout);
    
    return EXIT_SUCCESS;
}