CC=g++
CPPFLAGS=-c -g -std=c++0x -Wall -Wextra -pedantic
LDFLAGS = -lsfml-graphics -lsfml-window -lsfml-system
SOURCES=map.hpp map.cpp  main.cpp game.hpp game.cpp weapon.hpp weapon.cpp calc.hpp calc.cpp player.hpp player.cpp menu.cpp
OBJECTS=main.o game.o map.o player.o calc.o weapon.o menu.o
EXECUTABLE=Artillery

all: $(SOURCES) $(EXECUTABLE)
  
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

.cc.o:
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

run: Artillery
	./Artillery
