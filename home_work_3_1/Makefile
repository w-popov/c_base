# $@ - имя цели
# $< - имя первой зависимости
# $^ - имена всех зависимостей
# $* - имя цели без расширения
# %.o: %.c %.h - правило для компиляции .c в .o
# %: %.c — если нужно создать файл с именем some_name из файла some_name.c
# $(wildcard ...) — функция поиска файлов по шаблону

# цвета
GREEN=\033[0;32m
PURPLE=\033[0;35m
BLUE=\033[0;34m
RED=\033[0;31m
NC=\033[0m 

CC=gcc
CFLAGS=-Wall -Wextra
STANDART=-std=c99
LDFLAGS=-lm
SOURCES = $(wildcard *.c)
TARGETS = $(SOURCES:.c=)

all: $(TARGETS)

%: %.c
	@echo "$(BLUE)Компиляция: $(NC)"
	$(CC) $(STANDART) $(CFLAGS) $< $(LDFLAGS) -o $@
	@echo "$(PURPLE)Добавление прав: $(NC)"
	chmod +x $@
	@ls -l $@
	@echo "\n"


run: all
	@for exe in $(TARGETS); do \
		if [ -x $$exe ]; then \
			echo "$(GREEN) Запуск $$exe:$(NC)"; \
			./$$exe; \
			echo ""; \
			echo ""; \
		else \
			echo "$(RED)Ошибка: $$exe не является исполняемым файлом$(NC)"; \
			ls -l $$exe; \
		fi \
	done


clean:
	rm -f $(TARGETS)

.PHONY: all run clean