INC="./inc"
FLAGS=-I$(INC)
OMPFLAG=-fopenmp
CC=gcc

knn: knn.c
	$(CC) knn.c -o knn -lm $(OMPFLAG)

knn_parallel: knn_parallel.c
	$(CC) knn_parallel.c -o knn_parallel -lm $(OMPFLAG)

debug: knn.c
	$(CC) -g knn.c -o debug -lm $(OMPFLAG)

clean:
	rm main
