
SOURCES=zsniff.cpp zdata.cpp
CFLAGS=-std=c++0x
CC=g++

all:
	$(CC) $(CFLAGS) -o zsniff $(SOURCES) 


