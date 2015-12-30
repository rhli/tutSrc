CC = gcc
CFLAG = -shared -fPIC -I$(HOME)/include -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux

all: lib/libjnitest.so

lib/libjnitest.so : nativeSrc/JniTest.c 
	$(CC) -o $@ $(CFLAG) $<

clean:
	rm lib/libjnitest.so

