CC=gcc
CFLAGS=-c -Wall
LDFLAGS=
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=main
SOURCES=main.c add_num.c

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)