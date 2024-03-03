all: main.c ws
	gcc -Wall -Wextra -I./wsServer/include -c main.c
	gcc  main.o -L./wsServer -pthread -lws -lm -o gps

ws:
	make -C wsServer

clean:
	rm -rf *.o gps
	make -C wsServer clean