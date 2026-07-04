#include <stdio.h>
#include <stdlib.h>
#include "temp_api.h"
#include "csv_parse.h"

struct TemperatureStats stats_array[MAX_SIZE_ARRAY];

int main(void)
{
    const char* filename = "temperature_small.csv";
    struct TemperatureStats* parsed_stats = parse_csv(filename, stats_array, MAX_SIZE_ARRAY);
    if (parsed_stats == NULL)
    {
        fprintf(stderr, RED "Error parsing CSV file: %s\n" RESET, filename);
        return EXIT_FAILURE;
    }
    print_temperature_stats_array(parsed_stats, MAX_SIZE_ARRAY);

    return EXIT_SUCCESS;
}
