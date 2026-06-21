/**
 * ДЗ-9. Си базовый уровень. гр.Д01-134 Попов. В.Г
 * F16: Известно, что шахматная доска имеет 
 *      размерность 8х8 и состоит из клеток 2х цветов...
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/**
 * A4 => [col, row]
 * Если: row - чет и  | col - нечет, то white
 *                    | col - чет, то black
 * 
 * Если row - нечет и | col - нечет, то black
 *                    | col - чет, то white
 */
const char* black_or_white (char* coords)
{
    if (!coords || strlen(coords) != 2 || !(isalpha(coords[0]) && isdigit(coords[1])))
    {
        return "Error coordinates\n";
    }
    /* Буква маленькая */
    coords[0] = tolower(coords[0]);
    if ( !(coords[0] >= 'a' && coords[0] <= 'h') )
    {
        return "Error colums coordinates\n";
    }
    
    enum { BOARD = 9 };
    char cols[BOARD] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', '\0'};
    char *pt_col = strchr(cols, coords[0]);
    if (pt_col == NULL)
    {
        return "Error searching index\n";
    }
    /* Координаты col, row: */
    int col = pt_col - cols;
    int row = coords[1] - '0';
    if (row % 2)
    {
        if (col % 2)
        {
            return "WHITE";
        }
        return "BLACK";
    }
    else if (!(row % 2))
    {
        if (col % 2)
        {
            return "BLACK";
        }
    }
    return "WHITE";
}

#ifndef TEST_DEF_HW9
int main (void)
{
    enum { SIZE_ARR = 4 };
    char coordinates[SIZE_ARR] = {'\0'};
    if (scanf("%s", coordinates) != 1)
    {
        return EXIT_FAILURE;
    }
    printf("%s\n", black_or_white(coordinates));

    return EXIT_SUCCESS;
}
#endif