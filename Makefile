all: main.c
	gcc -Wall -Wextra -I./wsServer/include -c main.c
	make -C wsServer
	gcc  main.o -L./wsServer -pthread -lws -lm -o gps

clean:
	rm -rf *.o gps
	make -C wsServer clean