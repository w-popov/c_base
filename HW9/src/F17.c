/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F17: Составить функцию которая находит след матрицы 
 *      в двумерном массиве. Показать пример ее 
 *      работы на матрице из 5 на 5 элементов.
 */
#include <stdio.h>
#include <stdlib.h>

enum { SIZE_MATIX = 5 };

int trace_of_matrix (int matrix[][SIZE_MATIX], const int size)
{
    int trace_sum = 0;
    for (int i = 0; i < size; ++i)
    {
        trace_sum += matrix[i][i];
    }
    return trace_sum;
}

/**
 * Заполнить двумерный массив
 * @return Указатель на двумерный массив
 */
int (*fill_matrix_F17 (int matrix[][SIZE_MATIX], const int size))[SIZE_MATIX]
{
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

#ifndef TEST_DEF_HW9
int main (void)
{
    int array[SIZE_MATIX][SIZE_MATIX] = {0};
    int (*matrix)[SIZE_MATIX] = fill_matrix_F17(array, SIZE_MATIX);
    if (matrix == NULL)
    {
        printf("Fail fill matrix\n");
        return EXIT_FAILURE;
    }
    printf("%d\n", trace_of_matrix(matrix, SIZE_MATIX));
    return EXIT_SUCCESS;
}
#endif