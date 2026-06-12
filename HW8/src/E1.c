/**
 * ДЗ-8. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * E: 
 */
#include <stdio.h>
#include <stdlib.h>

int test_sum(void)
{
    int a = 0, b = 0;
    scanf("%d %d", &a, &b);
    return a + b;
}

#ifndef TEST_DEF_HW8
int main (void)
{
    int input_number = 0;
    scanf("%d\n", &input_number);
    
    return EXIT_SUCCESS;
}
#endif