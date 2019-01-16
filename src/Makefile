LIB = ../lib/
BIT_ARR_PATH = $(LIB)BitArray
BIT_ARR_ARGS = -I $(BIT_ARR_PATH) -L $(BIT_ARR_PATH) -lbitarr


all: codage packet

codage:	codage.c
	gcc -o codage.out codage.c random.c $(BIT_ARR_ARGS)

packet: packet.c
	gcc -c packet.c $(BIT_ARR_ARGS)

random: random.c
	gcc -c random.c $(BIT_ARR_ARGS)

libs:
	$(MAKE) -C $(LIB)

clean:
	rm *.o
	rm *.out

.PHONY: all clean libs