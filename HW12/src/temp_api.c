#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "temp_api.h"
#include "parcer_csv.h"


/* --------------------- вывод статистики по каждому месяцу: ------------*/ 
/**
 * Среднемесячная температура 
 */
int16_t average_monthly_temperature (struct TemperatureStats *tarr, uint16_t month)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument average_monthly_temperature()\n" RESET);
        return INT16_MIN;
    }
    if (month < 1 || month > 12)
    {
        perror(RED "Error! the month number must be between 1 and 12\n" RESET);
        return INT16_MIN;
    }
    long long int sum = 0;
    long long int count = 0;
    for (size_t i = 0; i < MAX_SIZE_ARRAY; ++i)
    {
        if ((uint16_t)tarr[i].month == month )
        {
            sum += (int16_t)tarr[i].temperature;
            count++;
        }
    }
    if (count == 0)
    {
        return INT16_MIN;
    }
    return (int16_t)(sum / count);
}

/**
 * Минимальная температура в текущем месяце 
 */
int16_t min_temperature_current_month (struct TemperatureStats *tarr, uint16_t month)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument min_temperature_current_month()\n" RESET);
        return INT16_MIN;
    }
    if (month < 1 || month > 12)
    {
        perror(RED "Error! the month number must be between 1 and 12\n" RESET);
        return INT16_MIN;
    }
    int16_t min_temp = INT16_MAX;
    for (size_t i = 0; i < MAX_SIZE_ARRAY; ++i)
    {
        if ((uint16_t)tarr[i].month == month)
        {
            min_temp = tarr[i].temperature;
        }
    }
    return min_temp;
}

/**
 *  Максимальная температура в текущем месяце 
 */
int16_t max_temperature_current_month (struct TemperatureStats *tarr, uint16_t month)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument max_temperature_current_month()\n" RESET);
        return INT16_MIN;
    }
    if (month < 1 || month > 12)
    {
        perror(RED "Error! the month number must be between 1 and 12\n" RESET);
        return INT16_MIN;
    }
    int16_t max_temp = INT16_MIN;
    for (size_t i = 0; i < MAX_SIZE_ARRAY; ++i)
    {
        if ((uint16_t)tarr[i].month == month)
        {
            max_temp = (int16_t)tarr[i].temperature;
        }
    }
    return max_temp;
}

/* --------------------- вывод статистику за год: ----------------------*/
/** 
 * Среднегодовая температура 
 */
int16_t average_annual_temperature (struct TemperatureStats *tarr, uint16_t year)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument average_annual_temperature()\n" RESET);
        return INT16_MIN;
    }
    if (year < 1950 || year > 2100)
    {
        perror(RED "Error! the year number must be between 1950 and 2100\n" RESET);
        return INT16_MIN;
    }
    long long int sum = 0;
    long long int count = 0;
    for (size_t i = 0; i < MAX_SIZE_ARRAY; ++i)
    {
        if ((uint16_t)tarr[i].year == year)
        {
            sum += (int16_t)tarr[i].temperature;
            count++;
        }
    }
    if (count == 0)
    {
        return INT16_MIN;
    }
    return (int16_t)(sum / count);
}

/**
 *  Минимальная температура 
 */
int16_t minimum_temperature (struct TemperatureStats *tarr, uint16_t year)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument minimum_temperature()\n" RESET);
        return INT16_MIN;
    }
    if (year < 1950 || year > 2100)
    {
        perror(RED "Error! the year number must be between 1950 and 2100\n" RESET);
        return INT16_MIN;
    }
    int16_t min_temp = INT16_MAX;
    for (size_t i = 0; i < MAX_SIZE_ARRAY; ++i)
    {
        if ((uint16_t)tarr[i].year == year)
        {
            min_temp = (int16_t)tarr[i].temperature;
        }
    }
    return min_temp;
}

/**
 *  Максимальная температура 
 */
int16_t maximum_temperature (struct TemperatureStats *tarr, uint16_t year)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument maximum_temperature()\n" RESET);
        return INT16_MIN;
    }
    if (year < 1950 || year > 2100)
    {
        perror(RED "Error! the year number must be between 1950 and 2100\n" RESET);
        return INT16_MIN;
    }
    int16_t max_temp = INT16_MIN;
    for (size_t i = 0; i < MAX_SIZE_ARRAY; ++i)
    {
        if ((uint16_t)tarr[i].year == year)
        {
            max_temp = (int16_t)tarr[i].temperature;
        }
    }
    return max_temp;
}

/**
 * Вывод массива температур
 */
void print_temperature_stats_array (struct TemperatureStats *tarr, size_t size)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument print_temperature_stats_array()\n" RESET);
        return;
    }

    setvbuf(stdout, NULL, _IOFBF, 65536);
    for (size_t i = 0; i < size; ++i)
    {
        printf("%6zu]  Year: %4u,  Month: %2u,  Day: %2u,  Hours: %2u,  Minutes: %2u,  Temperature: %2d\n",
            i + 1, tarr[i].year, tarr[i].month, tarr[i].day, tarr[i].hours, tarr[i].minutes, tarr[i].temperature);
    }

    fflush(stdout);
}

void show_help(void)
{
    printf(YELLOW"======================== help: ============================\n");
    printf("    -h Справка\n");
    printf("    -f Имя файла .csv\n");
    printf("    -m Вывод статистики по номеру месяца (если задан ключ)\n");
    printf("==============================================================\n\n"RESET);
}

/**
 * Экран открытия файла
 */
FILE* show_open_file_status (const char *file_name, int64_t *filesize)
{
    FILE *file = open_file(file_name, filesize);
    int64_t fsize = -1;

    #if defined(_WIN32) || defined(_WIN64)
    struct _stati64 statbuf;
    if (_stati64(filename, &statbuf) == 0)
    {
        fsize = statbuf.st_size;
    }
    #else 
    struct stat statbuf;
    if (stat(file_name, &statbuf) == 0)
    {
        fsize = statbuf.st_size;
    }
    #endif

    if (file == NULL || fsize == -1)
    {
        printf(RED_BOLD"Ошибка открытия файла: %s\n\n"RESET, file_name);
        return NULL;
    }

    *filesize = fsize;    
    double mb_size = (double)fsize / (1024 * 1024);
    printf("══════════════════ Открытие файла .csv: ═══════════════════\n");
    printf("\tфайл:   %s\n", file_name);
    printf("\tразмер: %.2f Mb\n", mb_size);
    printf("═══════════════════════════════════════════════════════════\n");
    return file;
}

