/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C12: Составить функцию, которая вычисляет
 *      синус как сумму ряда (с точностью 0.001)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

/**
 * Из градусов в радианы
 */
float to_radian (int deg)
{
    return ((float)deg * PI) / 180;
}

float sinus (float arg)
{
    enum { ONE = 1, TWO = 2 };
    const float PRECISION = 0.001f;
    float x = to_radian(arg);

    while (x > PI)
    {
        x -= (float)TWO * PI;
    }
    while (x < -PI)
    {
        x += (float)TWO * PI;
    }
    float sum = 0.0f, current = x;
    for (int i = ONE; fabs(current) >= PRECISION; i += TWO)
    {
        sum += current;
        current = -current * x * x / ((i + ONE) * (i + TWO));
    }
    return sum;
}

#ifndef TEST_DEF_HW6
int main (void)
{
    float input_num = 0.0;
    scanf("%f", &input_num);
    printf("%.3f\n", sinus(input_num));
    return EXIT_SUCCESS;
}
#endif