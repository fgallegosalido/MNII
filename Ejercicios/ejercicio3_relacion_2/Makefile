SHELL = /bin/bash

BIN = bin
SRC = source
INC = include

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic -I$(INC)

all: $(BIN)/ejercicio3

$(BIN)/ejercicio3: $(SRC)/ejercicio3.cpp $(INC)/polynomial.hpp $(INC)/aux.hpp
	@mkdir $(BIN) 2> /dev/null
	@echo "Compilando $@........"
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	-rm -r $(BIN) 2> /dev/null
