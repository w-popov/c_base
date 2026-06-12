##### Склонировать каталог с ДЗ (windows)
* Создать где-то каталог, зайти в него, открыть в нем терминал.
* Ввести в терминале: 
```Bash
git clone --no-checkout https://github.com/w-popov/c_base.git
```
* Далее:
```Bash
cd c_base/
git sparse-checkout set HW8
git checkout
cd HW8/
```

##### Задать разные компиляторы и разные сборщики:
Для windows с [Ninja](https://github.com/ninja-build/ninja/releases):
```Bash
# В корне каталога:
cmake -G "Ninja" -B build

cmake --build build

# --build build --verbose или 
# --build build -v вывод подробной информации в процессе сборки
```
Либо с указанием компиляторов:
```Bash
# Из корня проекта
cmake -B build \
    -G "Ninja" \
    -DCMAKE_C_COMPILER=D:/Qt/Tools/mingw1310_64/bin/gcc.exe \
    -DCMAKE_CXX_COMPILER=D:/Qt/Tools/mingw1310_64/bin/g++.exe
   
    # C:/mingw64/bin/gcc.exe — может быть другой путь к gcc и g++
```
Для Linux с "Ninja" (по умолчанию "Unix Makefiles" и gcc — что вполне достаточно):
```Bash
# Установить, если его нет
sudo apt-get update && sudo apt-get install -y cmake ninja-build build-essential
# Сборка из корня проекта
cmake -G "Ninja" -B build
cmake --build build
```
##### Сборка домашних заданий:
* Перейти в корневой каталог с CMakeLists.txt
* Создать каталог build (можно не создавать — сгенерируется автоматически) ```mkdir build```
1. Запуск конфигурации CMake: ```cmake -B build```
2. Запуск сборки:             ```cmake --build build``` (сборка в каталог build)
* Для запуска тестов перейти в build ```cd build``` и запустить тесты: ввести ```ctest --rerun-failed --output-on-failure``` enter, 
  либо так (находясь в корневом каталоге): ```./runtests/hwtests ```
* Для запуска отдельных файлов ДЗ перейти в в bin ```cd bin``` и запустить например ```./E1```
* Очистить build (обычно нужно полная чистка build после серьезных изменений в CMakeLists.txt):
    1. ```rm -rf build/*  ``` либо
    2. ```cmake --build build --target clean``` 
