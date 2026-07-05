#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <math.h>
#include "temp_api.h"
#include "csv_parse.h"
#ifdef _WIN32
#include <windows.h>
#endif

struct TemperatureStats stats_array[MAX_SIZE_ARRAY];
char filename[FILE_NAME_PATH] = {'\0'};

int main(void)
{
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    #else
    setlocale(LC_ALL, "en_US.UTF-8");
    #endif

    show_select_file_menu(filename);
    
    struct CsvParseStatus* parse_status = parse_csv(filename, stats_array, MAX_SIZE_ARRAY);
    if (parse_status == NULL)
    {
        return EXIT_FAILURE;
    }
    size_t error_count = errors_numbers(parse_status, MAX_SIZE_ARRAY);
    if (error_count > 0)
    {
        fprintf(stderr, RED_BOLD "Количество ошибок при парсинге: %zu\n" RESET, error_count);
        print_errors(parse_status, MAX_SIZE_ARRAY);
    }
    printf(GREEN "Парсинг CSV-файла выполнен\n" RESET);
    
    print_temperature_stats_array(stats_array, 12);

    int mode = -1;
    
    while (mode)
    {
        show_menu();
        if (scanf("%d", &mode) != 1)
        {
            break;
        }

        switch (mode)
        {
        case 0:
            break;
        case 1:
            {
                uint16_t month = 0;
                printf("Введите № месяца: ");
                scanf("%hu", &month);
                int16_t result = average_monthly_temperature(stats_array, month);
                print_result(result);
                break;
            }
        case 2:
            {
                uint16_t month = 0;
                printf("Введите № месяца: ");
                scanf("%hu", &month);
                int16_t result = min_temperature_current_month(stats_array, month);
                print_result(result);
                break;
            }
        case 3:
            {
                uint16_t month = 0;
                printf("Введите № месяца: ");
                scanf("%hu", &month);
                int16_t result = max_temperature_current_month(stats_array, month);
                print_result(result);
                break;
            }
        case 4:
            {
                uint16_t year = 0;
                printf("Введите год: ");
                scanf("%hu", &year);
                int16_t result = average_annual_temperature(stats_array, year);
                print_result(result);
                break;
            }
        case 5:
            {
                uint16_t year = 0;
                printf("Введите год: ");
                scanf("%hu", &year);
                int16_t result = minimum_temperature(stats_array, year);
                print_result(result);
                break;
            }
        case 6:
            {
                uint16_t year = 0;
                printf("Введите год: ");
                scanf("%hu", &year);
                int16_t result = maximum_temperature(stats_array, year);
                print_result(result);
                break;
            }
        
        default:
            break;
        }
    }
    
    free(parse_status);
    return EXIT_SUCCESS;
}
// temperature_big.csv