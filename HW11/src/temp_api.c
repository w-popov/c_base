#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "temp_api.h"



/* --------------------- вывод статистики по каждому месяцу: ------------*/ 
/**
 * Среднемесячная температура 
 */
float average_monthly_temperature (struct TemperatureStats *tarr, uint16_t month)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument average_monthly_temperature()\n" RESET);
        return NAN;
    }
    if (month < 1 || month > 12)
    {
        perror(RED "Error! the month number must be between 1 and 12\n" RESET);
        return NAN;
    }

    return 0.0;
}

/**
 * Минимальная температура в текущем месяце 
 */
float min_temperature_current_month (struct TemperatureStats *tarr, uint16_t month)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument min_temperature_current_month()\n" RESET);
        return NAN;
    }
    if (month < 1 || month > 12)
    {
        perror(RED "Error! the month number must be between 1 and 12\n" RESET);
        return NAN;
    }

    return 0.0;
}

/**
 *  Максимальная температура в текущем месяце 
 */
float max_temperature_current_month (struct TemperatureStats *tarr, uint16_t month)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument max_temperature_current_month()\n" RESET);
        return NAN;
    }
    if (month < 1 || month > 12)
    {
        perror(RED "Error! the month number must be between 1 and 12\n" RESET);
        return NAN;
    }

    return 0.0;
}

/* --------------------- вывод статистику за год: ----------------------*/
/** 
 * Среднегодовая температура 
 */
float average_annual_temperature (struct TemperatureStats *tarr, uint16_t year)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument average_annual_temperature()\n" RESET);
        return NAN;
    }
    if (year < 1950 || year > 2100)
    {
        perror(RED "Error! the year number must be between 1950 and 2100\n" RESET);
        return NAN;
    }

    return 0.0;
}

/**
 *  Минимальная температура 
 */
float minimum_temperature (struct TemperatureStats *tarr, uint16_t year)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument minimum_temperature()\n" RESET);
        return NAN;
    }
    if (year < 1950 || year > 2100)
    {
        perror(RED "Error! the year number must be between 1950 and 2100\n" RESET);
        return NAN;
    }

    return 0.0;
}

/**
 *  Максимальная температура 
 */
float maximum_temperature (struct TemperatureStats *tarr, uint16_t year)
{
    if (tarr == NULL)
    {
        perror(RED "Error! zero pointer argument maximum_temperature()\n" RESET);
        return NAN;
    }
    if (year < 1950 || year > 2100)
    {
        perror(RED "Error! the year number must be between 1950 and 2100\n" RESET);
        return NAN;
    }

    return 0.0;
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
        printf(GREEN "%zu] Year: %u, Month: %u, Day: %u, Hours: %u, Minutes: %u, Temperature: %.2f\n" RESET,
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
