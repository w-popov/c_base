/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F19: Определить количество положительных элементов 
 *      квадратной матрицы, превышающих по величине среднее 
 *      арифметическое всех элементов главной диагонали. 
 *      Реализовать функцию среднее арифметическое главной диагонали.
 */
#include <stdio.h>
#include <stdlib.h>

enum { SIZE_MATRIX = 5 };

/**
 * Заполнить матрицу
 * @return Указатель на 2-мерный массив
*/
int ( *fill_matrix_F19 (int matrix[][SIZE_MATRIX], const int size)) [SIZE_MATRIX]
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

/**
 * Массив указателей на элементы гл. диагонали
 */
int **main_diagonal (const int size_square_matrix)
{
    int **diagonal = calloc((size_t)size_square_matrix, sizeof(int*));
    if (diagonal == NULL)
    {
        perror("Fail allocation memory in fun main_diagonal()\n");
        return NULL;
    }
    return diagonal;
}

void delete_diagonal (int **pt)
{
    if (pt)
    {
        free(pt);
    }
}

/**
 * Среднее арифметическое гл. диагонали
 */
int mean_main_diagonal (int **diagonal, const int size_diagonal)
{
    if (diagonal == NULL)
    {
        perror("**diagonal pointer is NULL in fun mean_main_diagonal()\n");
        return 0;
    }
    int sum = 0, **pt_diagonal = diagonal;
    while (pt_diagonal < diagonal + size_diagonal)
    {
        sum += **pt_diagonal;
        ++pt_diagonal;
    }    
    return sum / size_diagonal;
}

/**
 * Количество положительных элементов матрицы превышающих по величине среднее 
 * арифметическое всех элементов главной диагонали
 */
int number_of_positive_elements (int matrix[][SIZE_MATRIX], const int size_matrix, int mean)
{
    int count = 0;
    for (int row = 0; row < size_matrix; ++row)
    {
        for (int col = 0; col < size_matrix; ++col)
        {
            count += matrix[row][col] > 0 && matrix[row][col] > mean ? 1 : 0;
        }
    }
    return count;
}

/**
 * Решение задания
 */
int above_average (int matrix[][SIZE_MATRIX], const int size_matrix)
{
    int **diagonal = main_diagonal(size_matrix);
    if (diagonal == NULL)
    {
        perror("**diagonal pointer is NULL in fun above_average()\n");
        return 0;
    }
    for (int i = 0; i < size_matrix; ++i)
    {
        diagonal[i] = &matrix[i][i];
    }
    int mean = mean_main_diagonal(diagonal, size_matrix);
    delete_diagonal(diagonal);

    int count = number_of_positive_elements(matrix, size_matrix, mean);

    return count;
}

#ifndef TEST_DEF_HW9
int main (void)
{
    int matrix[SIZE_MATRIX][SIZE_MATRIX] = {0};
    if (fill_matrix_F19(matrix, SIZE_MATRIX) == NULL)
    {
        return EXIT_FAILURE;
    }
    int result = above_average(matrix, SIZE_MATRIX);
    printf("%d\n", result);

    return EXIT_SUCCESS;
}
#endif