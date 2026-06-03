/**
 * ДЗ-6. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * C13: Составить функцию, которая вычисляет косинус
 *      как сумму ряда (с точностью 0.001)
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

float cosinus (float arg)
{
    enum { ZERO, ONE, TWO };
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
    float sum = 0.0f, current = 1.0f;
    for (int i = ZERO; fabs(current) >= PRECISION; i += TWO)
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
    printf("%.3f\n", cosinus(input_num));
    return EXIT_SUCCESS;
}
#endif