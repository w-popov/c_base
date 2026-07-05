#ifndef _TEMP_API_H_
#define _TEMP_API_H_

#ifdef __cplusplus
extern "C" {
    #endif

#include <stdint.h>

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

#define MAX_SIZE_ARRAY  1000000    // Максимальный размер массива температур
#define FILE_NAME_PATH  2000        // Макс. размер пути/имени файла

struct TemperatureStats
{
    uint32_t year       :11;
    uint32_t month      :4;
    uint32_t day        :6;
    uint32_t hours      :5;
    uint32_t minutes    :6;
    int16_t temperature;
};

/* --------------------- вывод статистики по каждому месяцу: */ 
/* Среднемесячная температура */
int16_t average_monthly_temperature (struct TemperatureStats*, uint16_t);

/* Минимальная температура в текущем месяце */
int16_t min_temperature_current_month (struct TemperatureStats*, uint16_t);

/* Максимальная температура в текущем месяце */
int16_t max_temperature_current_month (struct TemperatureStats*, uint16_t);

/* --------------------- вывод статистику за год: */
/* Среднегодовая температура */
int16_t average_annual_temperature (struct TemperatureStats*, uint16_t);

/* Минимальная температура */
int16_t minimum_temperature (struct TemperatureStats*, uint16_t);

/* Максимальная температура */
int16_t maximum_temperature (struct TemperatureStats*, uint16_t);

/* Вывод массива температур */
void print_temperature_stats_array (struct TemperatureStats*, size_t size);

/* Меню вывода статистики */
void show_menu (void);

/* Меню выбора файла */
void show_select_file_menu (char file_name_buff[]);

#ifdef __cplusplus
}
#endif

#endif