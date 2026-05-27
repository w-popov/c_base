* Перейти в корневой каталог с CMakeLists.txt
* Создать каталог build ```mkdir build```
1. Запуск конфигурации CMake: ```cmake -B build```
2. Запуск генерации для системы сборки ```cmake --build build```
* Для запуска тестов перейти в build ```cd build``` и запустить тесты: ввести ```ctest``` enter, либо так: ```./hwtests ```
* Для запуска отдельных файлов ДЗ перейти в в build ```cd build``` и запустить например ```./B1```
* Очистить build находясь в каталоге HW5:
    1. ```rm -rf build/*  ``` либо
    2. ```cmake --build build --target clean``` 
