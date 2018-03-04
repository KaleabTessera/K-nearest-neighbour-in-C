INC="./inc"
FLAGS=-I$(INC)
OMPFLAG=-fopenmp
CC=gcc

main: main.c
	$(CC) main.c -o main -lm

debug: main.c
	$(CC) -g main.c -o debug -lm

clean:
	rm main
