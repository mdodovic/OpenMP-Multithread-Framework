# C compiler
CC = gcc
CC_FLAGS = -std=gnu99 -fopenmp -O2
CC_DEFS = -D LIFE_VISUAL

gameoflife: gameoflife.c
	$(CC) $(CC_FLAGS) gameoflife.c -o gameoflife 

visual: gameoflife.c
	$(CC) $(CC_DEFS) $(CC_FLAGS) gameoflife.c -o gameoflife 

clean:
	rm -f gameoflife
