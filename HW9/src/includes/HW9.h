#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/* F1 */
void sort_array(int size, int a[]);

/* F2 */
void sort_even_odd(int n, int a[]);

/* F5 */
int find_max_array(int size, int a[]);

/* F6 */
int is_two_same(int size, int a[]);

/* F7 */
int compression (int a[], int b[], int N);

/* F8 */
int fill_array_F8 (int *array, const int size);

/* F9 */
void swap_negmax_last(int size, int a[]);

/* F11 */
struct MinTwo
{
    int min1;
    int min2;
};

struct MinTwo find_2_minimum (int *array, const int size);

/* F12 */
void change_max_min(int size, int a[]);

/* F13 */
int count_between(int from, int to, int size, int a[]);

#ifdef __cplusplus
}
#endif