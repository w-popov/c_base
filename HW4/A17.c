/**
 * ДЗ-4. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * A17: Ввести номер месяца и вывести название времени года.
 *      Время года на английском: winter, spring, summer, autumn
 */

#include <stdio.h>
#include <stdlib.h>


unsigned season(unsigned month)
{
    return (month % 12) / 3;
}

/**
 * TEST_DEF_HW определяется в Makefile для компиляции тестов в ./tests. 
 * Компиляция без main()
 */
#ifndef TEST_DEF_HW
int main(void)
{
    #define SEASONS 4

    const char* seasons[SEASONS] = {"winter", "spring", "summer", "autumn"}; 
    unsigned month = 0, it_season = 0;
    scanf("%u", &month);
    it_season = season(month);
    printf("%s\n", seasons[it_season]);

    return EXIT_SUCCESS;
}
#endif