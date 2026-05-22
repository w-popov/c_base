В кратце:
При данной структуре каталогов запустить `gentests.py` с аргументами: `gentests folder-name cpp-tests-file` т-е:

`./gentests ../HW4 all`, где `HW4` каталог с Д.З. уровнем выше, `all` - имя файла тестов без приставки и расширения.
В каталоге HW4 (например), если он существует, будет создан подкаталог `/tests`, в котором будут сгенерирован файлы: `test_all.cpp, Makefile`. В каталоге с файлом `gentests.py` должен находиться `catch.hpp`. `catch.hpp` будет скопирован в каталог `/tests` (например в `HW4/tests/`).
[Ссылка на catch.hpp](https://github.com/catchorg/Catch2/releases/download/v2.13.10/catch.hpp)
Для тестирования кода из .c исходников Д.З. нужно в `test_all.cpp` определить (функции например) как:
```c 
extern "C"
{
    struct EquationLine { float k; float b; };
    struct EquationLine equation_line(int, int, int, int);
}
```
а сам `.c` файл:
```c
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
```
После написания тестов ввести `make tests` находясь в каталоге `/tests`
Ввести `make clean` для очистки от сгенерированных make файлов