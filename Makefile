all:
	gcc -g -Wall -std=c99 -pthread -fopenmp Lab3IO.c lab3.c -o main 
	gcc datagen.c Lab3IO.c -o datagen
	gcc serialtester.c Lab3IO.c -o serialtester -lm

clean:
	rm -f main datagen serialtester