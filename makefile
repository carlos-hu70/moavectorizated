CFLAGS      = -Wall -O3 -march=native
LFLAGS      = -Wall -O3 
CC      = gcc

OBJ     = obj/heap.o obj/moastarv.o obj/graph.o obj/main_moastarv.o


all: moastarv

moastarv:  $(OBJ)
	$(CC) $(LFLAGS) -o moastarv $(OBJ)


obj/%.o: src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f obj/*.o moastarv
