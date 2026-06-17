/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F13: 
 */
#include <stdio.h>
#include <stdlib.h>



#ifndef TEST_DEF_HW9
int main (void)
{
    enum { SIZE_ARR = 10 };
    int array[SIZE_ARR] = {0};
    printf("%d\n", array[0]);

    return EXIT_SUCCESS;
}
#endif