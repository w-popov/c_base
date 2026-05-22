#!/usr/bin/env python3

"""г
МФТИ "Пуск" 2026. Гр. Д01-134 Попов В.Г. СИ базовый уровень.
Генерация шаблонов для тестов домашних заданий. Используется C++, catch2 
"""

import os
import sys
from enum import StrEnum

TEMPLATE_CPP = """
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <vector>


TEST_CASE( "vectors can be sized and resized", "[vector]" )
{
    // For each section, vector v is anew:

    std::vector<int> v( 5 );

    REQUIRE( v.size() == 5 );
    REQUIRE( v.capacity() >= 5 );

    SECTION( "resizing bigger changes size and capacity" )
    {
        v.resize( 10 );

        REQUIRE( v.size() == 10 );
        REQUIRE( v.capacity() >= 10 );
    }
    SECTION( "resizing smaller changes size but not capacity" )
    {
        v.resize( 0 );

        REQUIRE( v.size() == 0 );
        REQUIRE( v.capacity() >= 5 );
    }
    SECTION( "reserving bigger changes capacity but not size" )
    {
        v.reserve( 10 );

        REQUIRE( v.size() == 5 );
        REQUIRE( v.capacity() >= 10 );
    }
    SECTION( "reserving smaller does not change size or capacity" )
    {
        v.reserve( 0 );

        REQUIRE( v.size() == 5 );
        REQUIRE( v.capacity() >= 5 );
    }
}
"""

TEMPLATE_MAKEFILE = """

# -MMD   Генерирует файл зависимостей (.d файл) во время компиляции
#        Отслеживает только пользовательские заголовки (не системные)
#        Файл .d содержит список всех заголовочных файлов, от которых зависит .cpp файл
#        Не прерывает компиляцию (в отличие от -MM)
#
# -MP     Добавляет пустые (фиктивные) цели для каждого заголовочного файла
#         Защищает от ошибок при удалении заголовочных файлов
#         Без этого флага удаление .h файла вызовет ошибку make
# $@ - имя цели
# $< - имя первой зависимости
# $^ - имена всех зависимостей
# $* - имя цели без расширения
# %.o: %.c %.h - правило для компиляции .c в .o
# %: %.c — если нужно создать файл с именем some_name из файла some_name.c
# $(wildcard ...) — функция поиска файлов по шаблону


# цвета
GREEN=\\033[0;32m
PURPLE=\\033[0;35m
BLUE=\\033[0;34m
RED=\\033[0;31m
NC=\\033[0m 

CC=gcc
CXX=g++
CXXFLAGS=-Wall -Wextra -Wpedantic
CFLAGS=-Wall -Wextra -Wpedantic
CXXSTANDARD=-std=c++20
CSTANDARD=-std=c99
LDFLAGS=-lm
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
TARGETS = $(SOURCES:.cpp=)
C_SRC=..
C_FILES=$(wildcard $(C_SRC)/*.c)
C_OBJS=$(notdir $(C_FILES:.c=.o))

# Макроопределение TEST_DEF_HW экранирует main в исходниках си для тестов
TEST_DEF = -DTEST_DEF_HW

# Отключить встроенные правила
.SUFFIXES:

all: $(TARGETS)

$(C_OBJS): %.o: $(C_SRC)/%.c
	@echo "$(BLUE)Компиляция C: $< -> $@$(NC)"
	$(CC) $(CSTANDARD) $(TEST_DEF) $(CFLAGS)  -c $< -o $@

# Правило для создания исполняемого файла из объектного
%: %.o $(C_OBJS)
	@echo "$(PURPLE)Линковка $@:$(NC)"
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
	@echo "$(PURPLE)Добавление прав: $(NC)"
	chmod +x $@
	@ls -l $@
	@echo ""

# Правило для компиляции .cpp в .o
%.o: %.cpp
	@clear
	@pwd
	@echo ""
	@echo "$(BLUE)Компиляция $< в $@:$(NC)"
	$(CXX) $(CXXSTANDARD) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Правило для включения зависимостей (файлов .d, созданных с помощью -MMD)
-include $(OBJECTS:.o=.d)

tests: all
	@for exe in $(TARGETS); do \\
		if [ -x $$exe ]; then \\
			echo "$(GREEN)---> Запуск тестов $$exe:$(NC)"; \\
			./$$exe; \\
			echo ""; \\
		else \\
			echo "$(RED)Ошибка: $$exe не является исполняемым файлом$(NC)"; \\
			ls -l $$exe; \\
		fi \\
	done

clean:
	rm -f $(TARGETS) *.d *.o

.PHONY: all tests clean

"""

class TermColors(StrEnum):
    YE = '\033[33m'
    GR = '\033[32m'
    BL = '\033[36m'
    RD = '\033[31m'
    NC = '\033[0m'
    CS = '\033[2J' # Clear screen

class Startup:
    args: tuple[str, ...]
    __name_hw: str = ''
    __name_files: list[str]
    
    def __init__(self, args: list[str]) -> None:
        print(f"""\n{TermColors.CS}{TermColors.YE}Генерация шаблонов для тестов домашних заданий.{TermColors.NC}""")
        try:
            self.args = tuple(os.path.splitext(arg)[0] for arg in self.__check_args(*args))
            print(self.args)
        except ValueError as e:
            print(e)
            sys.exit(1)


    def __check_args(self, *args: str) -> tuple[str, ...]:
        if not args:
            raise ValueError(f"\n{TermColors.RD}Нет аргументов для генерации шаблонов тестов.{TermColors.NC}")
        if len(args) < 3:
            raise ValueError(f"\n{TermColors.RD}Недостаточно аргументов! Аргументы: <путь/Имя каталога с дз>, <имена файлов без расширения .c>, ...\n{TermColors.NC}")
        return args[1:]
    
    
    def __create_test_folder(self) -> None:
        path = os.path.join(os.getcwd(), self.__name_hw)
        if not os.path.exists(path):
            print(f"{TermColors.RD}{self.__name_hw} не существует{TermColors.NC}")
            sys.exit(1)
        path = os.path.join(path, "tests")
        if os.path.exists(path):
            print(f"{TermColors.RD}Каталог tests уже существует в {self.__name_hw}{TermColors.NC}")
            print(f"{TermColors.BL}Удалить его? (y/n){TermColors.NC}")
            choice = input().strip().lower()
            if choice == 'y':
                # Удалить существующий каталог tests и его содержимое
                for root, dirs, files in os.walk(path, topdown=False):
                    for name in files:
                        os.remove(os.path.join(root, name))
                    for name in dirs:
                        os.rmdir(os.path.join(root, name))
                os.rmdir(path)
                print(f"{TermColors.GR}Каталог tests удалён.{TermColors.NC}")
                sys.exit(0)
            else:
                print(f"{TermColors.GR}Операция отменена. Каталог tests не будет создан.{TermColors.NC}")
                sys.exit(0)
        os.mkdir(path)
        print(f"{TermColors.GR}Каталог tests успешно создан в {self.__name_hw}{TermColors.NC}")
        

    def __find_catch2(self) -> str:
        # Искать catch.hpp в текущей директории
        for file in os.listdir(os.getcwd()):
            if file == "catch.hpp":
                return os.path.join(os.getcwd(), file)
        return ''
    
    
    def __copy_catch2(self, dest_path: str) -> None:
        catch2_path = self.__find_catch2()
        if not catch2_path:
            print(f"{TermColors.RD}catch.hpp не найден в текущей директории! Пожалуйста, поместите catch.hpp в эту директорию и запустите скрипт снова.{TermColors.NC}")
            sys.exit(1)
        else:
            with open(catch2_path, 'r') as src, open(dest_path, 'w') as dst:
                dst.write(src.read())
            print(f"{TermColors.GR}catch.hpp успешно скопирован в {dest_path}{TermColors.NC}")


    def __create_test_files(self) -> None:
        for file in self.__name_files:
            path = os.path.join(os.getcwd(), self.__name_hw, "tests", f"test_{file}.cpp")
            with open(path, 'w') as f:
                f.write(TEMPLATE_CPP)
        makefile_path = os.path.join(os.getcwd(), self.__name_hw, "tests", "Makefile")
        with open(makefile_path, 'w') as f:
            f.write(TEMPLATE_MAKEFILE)
        print(f"{TermColors.GR}Шаблоны тестов успешно созданы для файлов: {', '.join(self.__name_files)} и Makefile{TermColors.NC}")

    
    def init(self) -> None:        
        self.__name_hw = self.args[0]
        self.__name_files = [*self.args[1:]]
        self.__create_test_folder()
        catch2_path = self.__find_catch2()
        if not catch2_path:
            print(f"{TermColors.RD}catch.hpp не найден в текущей директории! Пожалуйста, поместите catch.hpp в эту директорию и запустите скрипт снова.{TermColors.NC}")
            sys.exit(1)
        else:
            # Скопировать catch.hpp в папку tests
            dest_path = os.path.join(os.getcwd(), self.__name_hw, "tests", "catch.hpp")
            if not os.path.exists(dest_path):
                self.__copy_catch2(dest_path)
        self.__create_test_files()


def main() -> None:
    startup = Startup(sys.argv)
    startup.init()

if __name__ == "__main__":
    main()