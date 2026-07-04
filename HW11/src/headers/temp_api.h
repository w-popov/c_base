#ifndef _TEMP_API_H_
#define _TEMP_API_H_

#ifdef __cplusplus
extern "C" {
    #endif

#include <stdint.h>

// Цвета
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"
#define RESET "\033[0m"

#define MAX_SIZE_ARRAY  100

struct TemperatureStats
{
    uint16_t year;
    uint16_t month;
    uint16_t day;
    uint16_t hours;
    uint16_t minutes;
    float temperature;
};

/* Указатель на ф-ции температур */
typedef float (*StatsFunction)(struct TemperatureStats*, uint16_t);

/* --------------------- вывод статистики по каждому месяцу: */ 
/* Среднемесячная температура */
float average_monthly_temperature (struct TemperatureStats*, uint16_t);

/* Минимальная температура в текущем месяце */
float min_temperature_current_month (struct TemperatureStats*, uint16_t);

/* Максимальная температура в текущем месяце */
float max_temperature_current_month (struct TemperatureStats*, uint16_t);

/* --------------------- вывод статистику за год: */
/* Среднегодовая температура */
float average_annual_temperature (struct TemperatureStats*, uint16_t);

/* Минимальная температура */
float minimum_temperature (struct TemperatureStats*, uint16_t);

/* Максимальная температура */
float maximum_temperature (struct TemperatureStats*, uint16_t);

/* Вывод массива температур */
void print_temperature_stats_array (struct TemperatureStats*, size_t size);

#ifdef __cplusplus
}
#endif

#endif