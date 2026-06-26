##### Склонировать каталог с ДЗ (windows)
* Создать где-то каталог, зайти в него, открыть в нем терминал.
* Ввести в терминале: 
```Bash
git clone --no-checkout https://github.com/w-popov/c_base.git
```
* Далее:
```Bash
cd c_base/
git sparse-checkout set HW10
git checkout
cd HW9/
```

##### Задать разные компиляторы и разные сборщики:
Для windows с [Ninja](https://github.com/ninja-build/ninja/releases):
```Bash
# В корне каталога (режим Debug):
cmake -G "Ninja" -B build -DCMAKE_BUILD_TYPE=Debug

cmake --build build

# В корне каталога (режим Release):
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
cmake -G "Ninja" -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```
##### Сборка домашних заданий:
* Перейти в корневой каталог с CMakeLists.txt
* Создать каталог build (можно не создавать — сгенерируется автоматически) ```mkdir build```
1. Запуск конфигурации CMake: ```cmake -B build```
2. Запуск сборки:             ```cmake --build build``` (сборка в каталог build)
* Для запуска тестов перейти в build ```cd build``` и запустить тесты: ввести ```ctest --rerun-failed --output-on-failure``` enter, 
  либо так (находясь в корневом каталоге): ```./runtests/hwtests ```
* Для запуска отдельных файлов ДЗ перейти в в bin ```cd bin``` и запустить например ```./G1```
* Очистить build (обычно нужно полная чистка build после серьезных изменений в CMakeLists.txt):
    1. ```rm -rf build/*  ``` либо
    2. ```cmake --build build --target clean``` 

##### GUI (FLTK 1.4.5 с++):
В дополнение для себя я решаю некоторые задачи с gui библиотекой FLTK 1.4.5: ```HW10/gui```, исполняемые файлы: ```HW10/bin/gui```. Алгоритм решения - код си решенной задачи из файла ```HW10/src/common.c``` обернутый классом c++. Все собирается вместе. 

<img width="1816" height="1079" alt="hwgui" src="https://github.com/user-attachments/assets/737eafe9-48c1-4ea6-9a13-c4a4c25c2b6b" />

