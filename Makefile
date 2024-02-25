all: main.c
	gcc -Wall -Wextra -I/usr/local/include -c main.c
	gcc -L/usr/local/lib main.o -lgsl -lgslcblas -lm -o gps 