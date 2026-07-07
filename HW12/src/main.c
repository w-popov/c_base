#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <unistd.h>
#include "temp_api.h"
#include "parcer_csv.h"
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

    FILE * file = NULL;
    int64_t filesize = -1;
    int key = -1;

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

    while ((key = getopt(argc,argv,"f:m::")) != -1)
    {
        switch (key)
        {
        case 'f':
            file = show_open_file_status(optarg, &filesize);
            cntx.file_size = filesize;
            result = parse_csv(&cntx, file);
            if (result == NULL)
            {
                perror(RED_BOLD"Ошибка парсинга!\n"RESET);
                break;
            }
            show_errors(result->errors_parse, result->array->size);
            break;
        
        default:
            break;
        }
    }
    
    svector_free(&t_array);
    svector_free(&err_array);
    return EXIT_SUCCESS;
}
// temperature_big.csv