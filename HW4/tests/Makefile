

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
	@echo ""
	@echo "$(BLUE)Компиляция $< в $@:$(NC)"
	$(CXX) $(CXXSTANDARD) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Правило для включения зависимостей (файлов .d, созданных с помощью -MMD)
-include $(OBJECTS:.o=.d)

tests: all
	@for exe in $(TARGETS); do \
		if [ -x $$exe ]; then \
			echo "$(GREEN)---> Запуск тестов $$exe:$(NC)"; \
			./$$exe; \
			echo ""; \
		else \
			echo "$(RED)Ошибка: $$exe не является исполняемым файлом$(NC)"; \
			ls -l $$exe; \
		fi \
	done

clean:
	rm -f $(TARGETS) *.d *.o

.PHONY: all tests clean

