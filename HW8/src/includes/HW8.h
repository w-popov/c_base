#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/* E1 */
float average_array (float *array, const int size_array);
/* E3 - - - - - - - - - - - - - - - - - - - - - - - - - */
struct MinMax
{
    unsigned max_index;
    int max;
    unsigned min_index;
    int min;
};
struct MinMax min_max_from_array (int* array, const int size);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/


#ifdef __cplusplus
}
#endif