# # ДЗ-11. Си базовый уровень. гр.Д01-134 Попов. В.Г
# #
# # debug режим. команда:   make
# # realese режим. команда: make DEBUG=0

# Имена исполняемых файлов
TARGET = main
TEST_TARGET = test

# Исходные каталоги проекта
SRC_DIR = src
HEAD_DIR = $(SRC_DIR)/headers
BUILD_BASE_DIR = build

# Каталоги для тестов
TEST_DIR = tests
TEST_HEAD_DIR = $(TEST_DIR)/headers
LIBS_DIR = libs

CC = gcc
CXX = g++
CSTANDARD = -std=c99
CXXSTANDARD = -std=c++20

# Режим сборки по умолчанию (1 - Debug, 0 - Release)
DEBUG ?= 1

# Базовые флаги компилятора (-MMD -MP для автоматического отслеживания изменений в .h)
COMMON_FLAGS = $(CSTANDARD) -Wall -Wextra -Wpedantic -I$(HEAD_DIR) -MMD -MP

# Настройка путей и флагов в зависимости от выбранного режима
ifeq ($(DEBUG), 1)
    CFLAGS = $(COMMON_FLAGS) -g -O0 -DDEBUG
    CXXFLAGS = $(CXXSTANDARD) -Wall -Wextra -Wpedantic -g -O0 -I$(TEST_HEAD_DIR) -I$(HEAD_DIR) -DDEBUG
    LDFLAGS = -lm
    
    # Каталоги для Debug режима
    BUILD_DIR   = $(BUILD_BASE_DIR)/debug
    OBJ_DIR     = $(BUILD_DIR)/obj
    EXECUTE_DIR = $(BUILD_DIR)/bin
    LIB_OBJ_DIR = $(BUILD_DIR)/$(LIBS_DIR)
else
    CFLAGS = $(COMMON_FLAGS) -O3 -flto -DNDEBUG
    CXXFLAGS = $(CXXSTANDARD) -Wall -Wextra -Wpedantic -O3 -flto -I$(TEST_HEAD_DIR) -I$(HEAD_DIR) -DNDEBUG
    LDFLAGS = -flto=auto -static-libgcc -static-libstdc++ -static -s -lm
    
    # Каталоги для Release режима
    BUILD_DIR   = $(BUILD_BASE_DIR)/release
    OBJ_DIR     = $(BUILD_DIR)/obj
    EXECUTE_DIR = $(BUILD_DIR)/bin
    LIB_OBJ_DIR = $(BUILD_DIR)/$(LIBS_DIR)
endif

# Поиск исходников проекта (.c) и генерация путей для объектных файлов (.o)
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# Поиск всех файлов тестов (*.cpp), исключая модуль инициализации catch_module.cpp
TEST_SRCS = $(filter-out $(TEST_DIR)/catch_module.cpp, $(wildcard $(TEST_DIR)/*.cpp))
TEST_OBJS = $(patsubst $(TEST_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(TEST_SRCS))

# Автоматическая генерация зависимостей (.d) для всех объектных файлов
DEPS = $(OBJS:.o=.d) $(TEST_OBJS:.o=.d) $(CATCH_MODULE_OBJ:.o=.d)

# Путь к скомпилированной библиотеки Catch (catch_module.o в папку libs)
CATCH_MODULE_OBJ = $(LIB_OBJ_DIR)/catch_module.o

# Финальные пути к исполняемым файлам
EXEC_TARGET = $(EXECUTE_DIR)/$(TARGET)
EXEC_TEST_TARGET = $(EXECUTE_DIR)/$(TEST_TARGET)

# Главная цель — собирает основную программу и тесты
all: $(EXEC_TARGET) $(EXEC_TEST_TARGET)

# Линковка основного исполняемого файла проекта
$(EXEC_TARGET): $(OBJS) | $(EXECUTE_DIR)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC_TARGET) $(LDFLAGS)
	cp $(EXEC_TARGET) ./$(TARGET)

# Линковка исполняемого файла тестов
# Собирает воедино библиотеку из каталога libs и все .o тестов
$(EXEC_TEST_TARGET): $(CATCH_MODULE_OBJ) $(TEST_OBJS) | $(EXECUTE_DIR)
	$(CXX) $(CXXFLAGS) $(CATCH_MODULE_OBJ) $(TEST_OBJS) -o $(EXEC_TEST_TARGET) $(LDFLAGS)

# Компиляция файлов .c проекта в файлы .o внутри каталога obj
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Компиляция файлов .cpp тестов в файлы .o внутри каталога obj
$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Сборка Catch в виде статической библиотеки в каталог libs
$(CATCH_MODULE_OBJ): $(TEST_DIR)/catch_module.cpp $(TEST_HEAD_DIR)/catch.hpp | $(LIB_OBJ_DIR)
	@echo "Компиляция Catch в каталог libs..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Создание каталога для объектных файлов проекта и тестов
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Создание каталога для объектных файлов библиотек (libs)
$(LIB_OBJ_DIR):
	mkdir -p $(LIB_OBJ_DIR)

# Создание каталога для бинарных файлов
$(EXECUTE_DIR):
	mkdir -p $(EXECUTE_DIR)

# Полная очистка
clean:
	rm -rf $(BUILD_BASE_DIR)
	rm -f ./$(TARGET)

# Подключение сгенерированных файлов зависимостей
-include $(DEPS)

.PHONY: all clean
