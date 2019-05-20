LIB = ../lib/
BIT_ARR_PATH = $(LIB)BitArray
BIT_ARR_ARGS = -I $(BIT_ARR_PATH) -L $(BIT_ARR_PATH) -lbitarr
RPI_ARGS = -lwiringPi

LIST = list.c
RANDOM = random.c boxmuller.c
TURBO = turbocode.c $(RANDOM)
LDPC = ldpc.c
TESTS = tests.c $(TURBO)
RPI = rpi.c
BASIC = basic.c
DEMO = demo.c

MAIN = codage.c $(TURBO) $(RPI) $(LIST) $(LDPC) $(BASIC) $(DEMO)
MAIN_ARGS = $(BIT_ARR_ARGS) $(RPI_ARGS) -lm

all: codage libs

codage: $(MAIN)
	gcc -o codage $(MAIN) $(MAIN_ARGS)

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