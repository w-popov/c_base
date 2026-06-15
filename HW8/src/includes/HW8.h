#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/* E1 */
float average_array (float *array, const int size_array);
/* E3 */
struct MinMax
{
    unsigned max_index;
    int max;
    unsigned min_index;
    int min;
};
struct MinMax min_max_from_array (int* array, const int size);
/* E4 */
int sum_two_max_from_array (int *array, const int size_array);
/* E5 */
int sum_positive_array (int *array, const int size_array);
/* E7 */
int *half_reverse_array (int *array, const int size_array);
/* E8 */
int *third_reverse_array (int *array, const int size_array);
/* E9 */
int *rshift_array (int *array, int shift, const int size);
/* E11 */
typedef int (*Compare)(int, int);
int compate_last_digit (int a, int b);
int *sort_X_array (int *array, const int size_array, Compare cmp);
/* E13 */
int *round_number_array (int *array, int *result_aray, const int size_array);
/* E16 */
struct DictInt
{
    int key;
    int value;
    struct DictInt *next;
};
void add_node_dict (struct DictInt **head, int key, int value);
int find_max_value (struct DictInt *head);
void delete_dict (struct DictInt *head);
int more_often (int *array, const int size_array);
/* E20 */
int build_greatest_number (int *digits_array, int size_num);
int the_greatest_number (int number);


#ifdef __cplusplus
}
#endif