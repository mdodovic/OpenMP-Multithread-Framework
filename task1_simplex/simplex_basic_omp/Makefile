# C compiler
CC = gcc
CC_FLAGS = -fopenmp -O2

mandelbrot: 
	$(CC) $(CC_FLAGS) simplex.c -o simplex -lm 

clean:
	rm -f simplex
