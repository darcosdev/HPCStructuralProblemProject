CC=g++
CXXFLAGS=-O3 -std=c++11
LDLIBS=-llapack -lblas -lboost_program_options
TARGET=hpcmain

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

.PHONY: clean run


task1: $(TARGET)
	./$(TARGET) --length 10 --nr_elements 24 --area 0.012 --mInertia 1.044e-05 --yMod 210e09 
	python task1.py

clean:
	rm -f $(TARGET) *.o *.txt
