/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * 
 */

#include <stdio.h>
#include <stdlib.h>

int uadd(unsigned a, unsigned b)
{
    return a + b;
}

#ifndef TEST_DEF_HW5
int main(void)
{
    printf("\n======> BUILD WITH CMAKE SUCCESS: %u\n", uadd(54, 8));

    return EXIT_SUCCESS;
}
#endif