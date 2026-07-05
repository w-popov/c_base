#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "temp_api.h"

static int valid_month(uint16_t month)
{
    return month <= 12;
}

static int valid_year(uint16_t year)
{
    return year >= 1950 && year <= 2100;
}

static int valid_temperature(int16_t temperature)
{
    return temperature > -99 && temperature < 100;
}


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
        if (tarr[i].month == month )
        {
            if (valid_month(tarr[i].month) && valid_temperature(tarr[i].temperature))
            {
                sum += tarr[i].temperature;
                count++;
            }
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
        if (tarr[i].month == month)
        {
            if (valid_temperature(tarr[i].temperature) && tarr[i].temperature < min_temp)
            {
                min_temp = tarr[i].temperature;
            }
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
        if (tarr[i].month == month)
        {
            if (valid_temperature(tarr[i].temperature) && tarr[i].temperature > max_temp)
            {
                max_temp = tarr[i].temperature;
            }
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
        if (tarr[i].year == year)
        {
            if (valid_year(tarr[i].year) && valid_temperature(tarr[i].temperature))
            {
                sum += tarr[i].temperature;
                count++;
            }
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
        if (tarr[i].year == year)
        {
            if (valid_temperature(tarr[i].temperature) && tarr[i].temperature < min_temp)
            {
                min_temp = tarr[i].temperature;
            }
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
        if (tarr[i].year == year)
        {
            if (valid_temperature(tarr[i].temperature) && tarr[i].temperature > max_temp)
            {
                max_temp = tarr[i].temperature;
            }
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

    for (size_t i = 0; i < size; ++i)
    {
        if (tarr[i].year == 0 && tarr[i].month == 0 && tarr[i].day == 0)
        {
            continue;
        }
        printf("%6zu]  Year: %4u,  Month: %2u,  Day: %2u,  Hours: %2u,  Minutes: %2u,  Temperature: %2d\n",
               i + 1,
               tarr[i].year,
               tarr[i].month,
               tarr[i].day,
               tarr[i].hours,
               tarr[i].minutes,
               tarr[i].temperature
            );
    }
}

void show_menu(void)
{
    printf(CYAN"======================== Введите: ============================\n");
    printf(" 0 - Выход\n");
    printf(" 1 - Среднемесячная температура\n");
    printf(" 2 - Минимальная температура в текущем месяце\n");
    printf(" 3 - Максимальная температура в текущем месяце\n");
    printf(" 4 - Среднегодовая температура\n");
    printf(" 5 - Минимальная годовая температура\n");
    printf(" 6 - Максимальная годовая температура\n");
    printf("==============================================================\n"RESET);

}

/* Меню выбора файла */
void show_select_file_menu (char file_name_buff[])
{
    printf("================== Введите имя файла .csv: ===================\n"MAGENTA);
    scanf("%[^\n]", file_name_buff);
    printf(RESET"Файл: %s\n", file_name_buff);
    printf("==============================================================\n\n");

}
