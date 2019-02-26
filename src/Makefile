LIB = ../lib/
BIT_ARR_PATH = $(LIB)BitArray
BIT_ARR_ARGS = -I $(BIT_ARR_PATH) -L $(BIT_ARR_PATH) -lbitarr

RANDOM = random.c boxmuller.c
TURBO = turbocode.c $(RANDOM)
TESTS = tests.c $(TURBO)

all: codage libs

codage:	codage.c random.c turbocode.c
	gcc -o codage codage.c $(TURBO) $(BIT_ARR_ARGS) -lm

packet: packet.c
	gcc -c packet.c $(BIT_ARR_ARGS)

random: random.c
	gcc -c random.c $(BIT_ARR_ARGS)

libs:
	@ cd $(LIB) && make

tests: $(TESTS)
	gcc -o tests $(TESTS) $(BIT_ARR_ARGS) -lm

clean:
	@ rm *.o
	@ echo Remove object files
	@ rm *.out
	@ echo Remove compiled files

.PHONY: all clean libs