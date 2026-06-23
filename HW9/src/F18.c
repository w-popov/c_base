/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F18: Дана целочисленная квадратная матрица 10 х 10. 
 *      реализовать алгоритм вычисления суммы максимальных 
 *      элементов из каждой строки. Напечатать значение этой суммы.
 *      Предполагается, что в каждой строке такой элемент единственный.
 *      Реализовать функцию поиска максимума в строке из 10 элементов
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int **new_matrix (const int size)
{
    if (size <= 0)
    {
        printf("Invalid size of matrix \n");
        return NULL;
    }
    int **matrix = calloc((size_t)size, sizeof(int*));
    if (matrix == NULL)
    {
        printf("Pointer is NULL \n");
        return NULL;
    }
    for (int i = 0; i < size; ++i)
    {
        matrix[i] = calloc((size_t)size, sizeof(int));
        if (matrix[i] == NULL)
        {
            printf("Fail allocation memory \n");
            for (int j = 0; j < i; ++j)
            {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

int free_matrix (int **matrix, const int size)
{
    if ( (matrix == NULL) || (size <= 0) )
    {
        printf("Invalid size or pointer is NULL\n");
        return 0;
    }
    for (int i = 0; i < size; ++i)
    {
        free(matrix[i]);
    }
    free(matrix);
    return 1;
}

int **fill_matrix_F18 (int **matrix, const int size)
{
    if ( (matrix == NULL) || (size <= 0) )
    {
        printf("Invalid size or pointer is NULL\n");
        return NULL;
    }
    for (int row = 0; row < size; ++row)
    {
        for (int col = 0; col < size; ++col)
        {
            if (scanf("%d", &matrix[row][col]) != 1)
            {
                return NULL;
            }
        }
    }
    return matrix;    
}

int get_max_value_row (int *row, const int size)
{
    if (size <= 0)
    {
        printf("Invalid size of matrix \n");
        return 0;
    }
    int max = INT32_MIN;
    for (int i = 0; i < size; ++i)
    {
        max = row[i] > max ? row[i] : max;
    }
    return max;
}

int sum_of_maxima (int **matrix, const int size)
{
    if ( (matrix == NULL) || (size <= 0) )
    {
        printf("Invalid size or pointer is NULL\n");
        return 0;
    }
    int sum = 0;
    for (int row = 0; row < size; ++row)
    {
        sum += get_max_value_row(matrix[row], size);
    }
    return sum;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { MATRIX_SIZE = 10 };
    
    int **matrix_10x10 = new_matrix(MATRIX_SIZE);
    if (matrix_10x10 == NULL)
    {
        return EXIT_FAILURE;
    }
    
    if (fill_matrix_F18(matrix_10x10, MATRIX_SIZE) == NULL)
    {
        printf("Input error\n");
        free_matrix(matrix_10x10, MATRIX_SIZE);
        return EXIT_FAILURE;
    }
    
    int sum = sum_of_maxima(matrix_10x10, MATRIX_SIZE);
    printf("%d\n", sum);
    
    free_matrix(matrix_10x10, MATRIX_SIZE);
    return EXIT_SUCCESS;
}
#endif