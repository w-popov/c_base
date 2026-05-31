/**
 * ДЗ-5. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * B17: Ввести натурально число и напечатать все числа 
 *      от 10 до введенного числа - у которых сумма цифр 
 *      равна произведению цифр
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW5.h"

typedef unsigned uint;

#ifndef TEST_DEF_HW5
int main(void)
{
    uint number = 0;
    scanf("%u", &number);
    UNIQE_PTR_U ptr_happy = happy_numbers(number);
    for (uint i = 0; i <= number; ++i)
    {
        if (ptr_happy.u_ptr[i])
        printf("%u ", ptr_happy.u_ptr[i]);
        else 
        {
            printf("\n");
            break;
        }
    }

    return EXIT_SUCCESS;
}
#endif


/**
 * Автоматическое освобождение аллоц. памяти
 * #define UNIQE_PTR_U __attribute__((cleanup(auto_free))) struct UniquePtr_u 
 */
void auto_free (struct UniquePtr_u* this)
{
    if (this->u_ptr) 
        free(this->u_ptr);
}

struct UniquePtr_u happy_numbers (uint number)
{
    enum { DIVIDER=10 };
    struct UniquePtr_u happy_arrays = { .u_ptr = (uint*)calloc( number + 1, sizeof(uint) ) };
    if (!happy_arrays.u_ptr)
    {
        printf("\033[33m Ошибка выделения памяти!(Allocation error.) \033[0m \n");
        exit(EXIT_FAILURE);
    }
    for (uint i = 0, index = 0, target = 10; i <= number; ++i, ++target)
    {
        uint summ = 0, mult = 1, _target = target;
        while (_target)
        {
            uint rest = _target % DIVIDER;
            summ += rest;
            mult *= rest;
            _target /= DIVIDER;
        }
        if (summ == mult)
            happy_arrays.u_ptr[index++] = target; 
    }
    return happy_arrays;
}