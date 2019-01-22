LIB = ../lib/
BIT_ARR_PATH = $(LIB)BitArray
BIT_ARR_ARGS = -I $(BIT_ARR_PATH) -L $(BIT_ARR_PATH) -lbitarr


all: codage libs

codage:	codage.c random.c turbocode.c
	gcc -o codage codage.c random.c turbocode.c $(BIT_ARR_ARGS)

packet: packet.c
	gcc -c packet.c $(BIT_ARR_ARGS)

random: random.c
	gcc -c random.c $(BIT_ARR_ARGS)

libs:
	@ cd $(LIB) && make

clean:
	@ rm *.o
	@ echo Remove object files
	@ rm *.out
	@ echo Remove compiled files

.PHONY: all clean libs