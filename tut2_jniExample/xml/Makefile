CC = g++

all: test

xmltest: src/test.cc src/tinyxml2.cpp src/tinyxml2.h
	$(CC) -o $@ $^

test: clean xmltest
	./xmltest
clean:
	rm -f *.o xmltest
