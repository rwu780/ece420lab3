all: datagen serialtester
	gcc -g -Wall -std=c99 -fopenmp Lab3IO.c lab3.c -o main 

baseline: datagen serialtester
	gcc -g -Wall -std=c99 -fopenmp Lab3IO.c Lab3_baseline_ompsolution.c -o main

datagen:
	gcc datagen.c Lab3IO.c -o datagen

serialtester:
	gcc serialtester.c Lab3IO.c -o serialtester -lm

clean:
	rm -f main datagen serialtester
