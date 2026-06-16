/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E17: Дан массив из 10 элементов. В массиве найти элементы, 
 *      которые в нем встречаются только один раз, и вывести их на экран. 
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
        if (scanf("%d", &array[i]) != 1)
        {
            printf("Error scanf\n");
            exit(EXIT_FAILURE);
        }
    }
    return array;
}

/**
 * Развернуть печать через рекурсию
 */
void print_just_once (struct DictInt *dict)
{
    if (dict == NULL)
    {
        return;
    }
    print_just_once(dict->next);
    if (!dict->value)
    {
        printf("%d ", dict->key);
    }
}

/**
 * Решение задания
*/
void just_once (int *array, const int size_array)
{
    struct DictInt *dict = NULL;
    for (int i = 0; i < size_array; ++i)
    {
        add_node_dict(&dict, array[i], 1);
    }
    print_just_once(dict);
    printf("\n");
    delete_dict(dict);
}

#ifndef TEST_DEF_HW8
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    just_once (
            fill_array(array, SIZE_ARR),
            SIZE_ARR
        );
    
    return EXIT_SUCCESS;
}
#endif
