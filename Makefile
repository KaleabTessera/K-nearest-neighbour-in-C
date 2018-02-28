INC="./inc"
FLAGS=-I$(INC)
OMPFLAG=-fopenmp
CC=gcc

main: main.c
	$(CC) main.c -o main -lm

clean:
	rm main
