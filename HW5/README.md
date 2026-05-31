##### Задать разные компиляторы и разные сборщики:
Для windows с [Ninja](https://github.com/ninja-build/ninja/releases):
```Bash
# Из корня проекта
cmake -B build \
    -G "Ninja" \
    -DCMAKE_C_COMPILER=C:/mingw64/bin/gcc.exe \
    -DCMAKE_CXX_COMPILER=C:/mingw64/bin/g++.exe
   
    # Путь к исходникам — текущая директория (.)
    # C:/mingw64/bin/gcc.exe -- может быть другой путь к gcc и g++
```
Для Linux с "Ninja" (по умолчанию "Unix Makefiles" и gcc):
```Bash
# Установить, если его нет
sudo apt-get update && sudo apt-get install -y cmake ninja-build build-essential
# Сборка из корня проекта
cmake -B build \
    -G "Ninja" \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++
```
##### Сборка домашних заданий:
* Перейти в корневой каталог с CMakeLists.txt
* Создать каталог build       ```mkdir build```
1. Запуск конфигурации CMake: ```cmake -B build```
2. Запуск сборки:             ```cmake --build build```
* Для запуска тестов перейти в build ```cd build``` и запустить тесты: ввести ```ctest --rerun-failed --output-on-failure``` enter, 
  либо так (находясь в корневом каталоге): ```./runtests/hwtests ```
* Для запуска отдельных файлов ДЗ перейти в в bin ```cd bin``` и запустить например ```./B1```
* Очистить build находясь в каталоге HW5:
    1. ```rm -rf build/*  ``` либо
    2. ```cmake --build build --target clean``` 
