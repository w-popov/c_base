# $@ - имя цели
# $< - имя первой зависимости
# $^ - имена всех зависимостей
# $* - имя цели без расширения
# %.o: %.cpp %.hpp - правило для компиляции .cpp в .o

CC=g++
STANDARD=-std=c++17
CFLAGS=-c -Wall
LDFLAGS=
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main
RUN=./main
CLEANOBJ=rm -rf *.o
SOURCES=main.cpp

all: $(SOURCES) $(EXECUTABLE)
	$(RUN) 
	$(CLEANOBJ)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(STANDARD) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(STANDARD) $(CFLAGS) $< -o $@

run: $(EXECUTABLE)
	$(RUN)

clean:
	$(CLEANOBJ) $(EXECUTABLE)