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

int main(void)
{
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    #else
    setlocale(LC_ALL, "en_US.UTF-8");
    #endif

    const char* filename = "temperature_big.csv";
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
    
    printf(GREEN "Парсинг CSV-файла прошел успешно\n" RESET);
    
    int mode = -1;
    
    while (mode)
    {
        show_menu();
        if (scanf("%d", &mode) != 1)
        {
            break;
        }

        switch ((int)mode)
        {
        case 0:
            break;
        case 1:
            {
                uint16_t month = 0;
                printf("Введите № месяца: ");
                scanf("%hu", &month);
                float result = average_monthly_temperature(stats_array, month);
                if (!isnan(result))
                {
                    printf(GREEN_BOLD"Результат: %d\n"RESET, (int)result);
                }
                break;
            }
        case 2:
            {
                uint16_t month = 0;
                printf("Введите № месяца: ");
                scanf("%hu", &month);
                float result = min_temperature_current_month(stats_array, month);
                if (!isnan(result))
                {
                    printf(GREEN_BOLD"Результат: %d\n"RESET, (int)result);
                }
                break;
            }
        case 3:
            {
                uint16_t month = 0;
                printf("Введите № месяца: ");
                scanf("%hu", &month);
                float result = max_temperature_current_month(stats_array, month);
                if (!isnan(result))
                {
                    printf(GREEN_BOLD"Результат: %d\n"RESET, (int)result);
                }
                break;
            }
        case 4:
            {
                uint16_t year = 0;
                printf("Введите год: ");
                scanf("%hu", &year);
                float result = average_annual_temperature(stats_array, year);
                if (!isnan(result))
                {
                    printf(GREEN_BOLD"Результат: %d\n"RESET, (int)result);
                }
                break;
            }
        case 5:
            {
                uint16_t year = 0;
                printf("Введите год: ");
                scanf("%hu", &year);
                float result = minimum_temperature(stats_array, year);
                if (!isnan(result))
                {
                    printf(GREEN_BOLD"Результат: %d\n"RESET, (int)result);
                }
                break;
            }
        case 6:
            {
                uint16_t year = 0;
                printf("Введите год: ");
                scanf("%hu", &year);
                float result = maximum_temperature(stats_array, year);
                if (!isnan(result))
                {
                    printf(GREEN_BOLD"Результат: %d\n"RESET, (int)result);
                }
                break;
            }
        
        default:
            break;
        }
    }
    

    free(parse_status);
    return EXIT_SUCCESS;
}
