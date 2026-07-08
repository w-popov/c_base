#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <unistd.h>
#include "temp_api.h"
#include "parser_csv.h"
#include "validate.h"

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#include <getopt.h>
#else
#include <unistd.h>
#endif

extern char *optarg; 
extern int optind, opterr, optopt; 

int main(int argc, char *argv[])
{
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    #else
    setlocale(LC_ALL, "en_US.UTF-8");
    #endif

    if (argc == 1)
    {
        show_help();
        return EXIT_SUCCESS;
    }

    FILE * file = NULL;
    int64_t filesize = -1;
    int key = -1;

    // Переменные для сохранения аргументов командной строки
    char *filename_arg = NULL;
    uint16_t month_arg = 0;     // 0 означает вывод за весь год
    int is_print = 0;           // вывод на экран?

    struct SVector t_array, err_array;
    struct SVector *array = svector_init(&t_array, sizeof(struct TemperatureStats), 0);
    struct SVector *errors_array = svector_init(&err_array, sizeof(struct ErrorParse), 0);
    CallbackProgressBar progressbar = print_progress_bar;
    CallbackWriteToArray write = write_to_array;
    struct ContextParser *result = NULL;

    struct ContextParser cntx = {
        .array = array,
        .delimiter = ";",
        .file_size = filesize,
        .errors_parse = errors_array,
        .length_field = 0,
        .current_column = 0,
        .current_row = 0,
        .nums_field = 5,
        .clb_progress = progressbar,
        .clb_write_to_arr = write
    };

    // разбор флагов и сохранение в параметры
    opterr = 0; 
    while ((key = getopt(argc, argv, "f:m:hp")) != -1)
    {
        switch (key)
        {
        case 'f':
            filename_arg = optarg;
            break;
        case 'p':
            
            is_print = 1;
            break;
        case 'm':
            
            month_arg = (uint16_t)atoi(optarg);
            if (month_arg < 1 || month_arg > 12)
            {
                fprintf(stderr, RED_BOLD "Ошибка: Неверный номер месяца '%s'. Ожидается от 1 до 12.\n" RESET, optarg);
                svector_free(&t_array);
                svector_free(&err_array);
                return EXIT_FAILURE;
            }
            break;
        case 'h':
            show_help();
            svector_free(&t_array);
            svector_free(&err_array);
            return EXIT_SUCCESS;
        case '?':
            if (optopt == 'f' || optopt == 'm')
            {
                fprintf(stderr, RED_BOLD "Ошибка: Флаг -%c требует указания значения.\n" RESET, optopt);
            }
            else
            {
                fprintf(stderr, RED_BOLD "Ошибка: Неизвестный флаг -%c.\n" RESET, optopt);
            }
            show_help();
            svector_free(&t_array);
            svector_free(&err_array);
            return EXIT_FAILURE;
        default:
            break;
        }
    }

    // Проверка наличия обязательного аргумента
    if (filename_arg == NULL)
    {
        fprintf(stderr, RED_BOLD "Ошибка: Не указан обязательный флаг -f (путь к CSV-файлу).\n" RESET);
        show_help();
        svector_free(&t_array);
        svector_free(&err_array);
        return EXIT_FAILURE;
    }

    // Открытие файла и выполнение парсинга
    file = show_open_file_status(filename_arg, &filesize);
    if (file == NULL)
    {
        svector_free(&t_array);
        svector_free(&err_array);
        return EXIT_FAILURE;
    }
    
    cntx.file_size = filesize;

    // парсить из файла
    struct ParseSource file_source = {
        .stream = file,
        .get_char = get_char_from_file,
        .get_pos = get_pos_from_file
    };

    result = parse_csv(&cntx, &file_source);
    fclose(file);
    file = NULL;

    if (result == NULL)
    {
        fprintf(stderr, RED_BOLD "Ошибка парсинга!\n" RESET);
        svector_free(&t_array);
        svector_free(&err_array);
        return EXIT_FAILURE;
    }

    show_errors(result->errors_parse, result->array->size);
    if (!is_print)
    {
        show_statistics((struct TemperatureStats*)array->data, array->size, month_arg);
    }
    else
    {
        print_temperature_stats_array((struct TemperatureStats*)array->data, array->size);
    }
    
    // Освобождение памяти
    svector_free(&t_array);
    svector_free(&err_array);

    #if defined(_WIN32) || defined(_WIN64)
    printf("\nНажмите Enter для выхода...\n");
    system("pause");
    #endif 

    return EXIT_SUCCESS;
}

// temperature_big.csv
// valgrind --leak-check=full --track-origins=yes ./main -f temperature_small.csv