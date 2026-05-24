/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A15: Определить уравнение прямой по координатам двух точек. Уравнение вида
 *      y=k*x+b. Ответ вида: два числа K,B в формате "%.2f %.2f"
 */

#include <stdio.h>
#include <stdlib.h>

struct EquationLine
{
    float k;
    float b;
};

struct EquationLine equation_line(int x1, int y1, int x2, int y2)
{
    float k = (float)(y2 - y1) / (x2 - x1);
    float b = (float)(y1) - (k * x1);
    return (struct EquationLine){k, b};
}


/**
 * TEST_DEF_HW определяется в Makefile для компиляции тестов в ./tests. 
 * Компиляция без main()
 */
#ifndef TEST_DEF_HW
int main(void)
{
    int x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    scanf("%d %d %d %d", &x1, &y1, &x2, &y2);
    struct EquationLine eq = equation_line(x1, y1, x2, y2);
    printf("%.2f %.2f\n", eq.k, eq.b);

    return EXIT_SUCCESS;
}
#endif