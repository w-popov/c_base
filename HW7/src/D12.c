/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D12: Дана монотонная последовательность, в которой каждое
 *      натуральное число k встречается ровно k раз: 1, 2, 2, 3, 3, 3,…
 *      Выведите первые n членов этой последовательности
 */
#include <stdio.h>
#include <stdlib.h>

/**
 * num - целевое кол-во чисел
 * it_num - текущее число
 * repeat - кол-во повторений
 */
void print_subsequence_rec (int num, int it_num, int repeat)
{
    if (!num)
    {
        return;
    }
    printf("%d ", it_num);

    if (repeat > 1)
    {
        print_subsequence_rec(num - 1, it_num, repeat - 1);
    }
    else
    {
        print_subsequence_rec(num - 1, it_num + 1, it_num + 1);
    }
}

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    print_subsequence_rec(input_number, 1, 1);
    printf("\n");
    return EXIT_SUCCESS;
}
#endif