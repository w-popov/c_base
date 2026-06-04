/**
 * ДЗ-7. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * D20: 
 */
#include <stdio.h>
#include <stdlib.h>
#include "HW7.h"

#ifndef TEST_DEF_HW7
int main (void)
{
    int input_number = 0;
    scanf("%d", &input_number);
    printf("%u\n", input_number);
    return EXIT_SUCCESS;
}
#endif