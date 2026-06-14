/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E16: Дан массив из 10 элементов. Определить, какое число 
 *      в массиве встречается чаще всего. Гарантируется, что такое число ровно 1. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/**
 * Упрощенный односвязный список - словарь.
 */
struct DictInt
{
    int key;
    int value;
    struct DictInt *next;
};

/**
 * Добавить узел. Если ключ сущесвует в списке, то новый узел
 * не добавляется, а к существующему узлу к полю value 
 * добавляется int value из параметра.
 */
void add_node_dict (struct DictInt **head, int key, int value)
{
    struct DictInt *find_node = *head;
    
    while (find_node != NULL)
    {
        if (find_node->key == key)
        {
            find_node->value += value;
            return;
        }
        find_node = find_node->next;
    }
    
    struct DictInt *node = malloc(sizeof(struct DictInt));
    node->key = key;
    node->value = 0;
    node->next = *head; 
    *head = node;
}

/**
 * Поиск чаще встречающегося значения
 */
int find_max_value (struct DictInt *head)
{
    int maxv = INT16_MIN, maxk = INT16_MIN;
    struct DictInt *current_dict = head;
    if (current_dict == NULL) {
        return 0;
    }
    while (current_dict != NULL)
    {
        if (current_dict->value > maxv)
        {
            maxv = current_dict->value;
            maxk = current_dict->key;
        }
        current_dict = current_dict->next;
    }
    return maxk;
}

/**
 * Удалить словарь. Освободить память.
 */
void delete_dict (struct DictInt *head)
{
    struct DictInt *bucket_dict = head;
    while (bucket_dict != NULL)
    {
        struct DictInt *next_node = bucket_dict->next;
        free(bucket_dict);
        bucket_dict = next_node;
    }
}
/* - - - - - - - - - - - - - - - - - - - - - - - -*/

int *fill_array (int *array, const int size_array)
{
    for (int i = 0; i < size_array; ++i)
    {
        scanf("%d", &array[i]);
    }
    return array;
}

/**
 * Решение задания
*/
int more_often (int *array, const int size_array)
{
    int often_num = 0;
    struct DictInt *dict = NULL;
    for (int i = 0; i < size_array; ++i)
    {
        add_node_dict(&dict, array[i], 1);
    }
    often_num = find_max_value(dict);
    delete_dict(dict);
    return often_num;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    int often = more_often (
                fill_array(array, SIZE_ARR), 
                SIZE_ARR
            );
    printf("%d\n", often);
    
    return EXIT_SUCCESS;
}
#endif