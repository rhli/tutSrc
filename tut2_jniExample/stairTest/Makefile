CC=gcc
CFLAG=-lgf_complete -I/usr/local/include
LDFLAG=-L/usr/local/lib

all:main

main:src/test.c stair.o NativeStair.o
	$(CC) -o $@ $^ $(CFLAG)

stair.o:src/stair.c src/stair.h 
	$(CC) -o $@ -c $< #$(CFLAG)

NativeStair.o:src/NativeStair.c src/NativeStair.h stair.o
	$(CC) -o $@ -c $< #$(CFLAG)

clean:
	rm *.o main
