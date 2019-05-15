LIB = ../lib/
BIT_ARR_PATH = $(LIB)BitArray
BIT_ARR_ARGS = -I $(BIT_ARR_PATH) -L $(BIT_ARR_PATH) -lbitarr
BMP_PATH = $(LIB)libbmp
BMP_ARGS = -I $(BMP_PATH) -L $(BMP_PATH) -lbmp
RPI_ARGS = -lwiringPi

LIST = list.c
RANDOM = random.c boxmuller.c
TURBO = turbocode.c $(RANDOM)
LDPC = ldpc.c
TESTS = tests.c $(TURBO)
RPI = rpi.c
BASIC = basic.c

MAIN = $(TURBO) $(RPI) $(LIST) $(LDPC) $(BASIC)
MAIN_ARGS = $(BIT_ARR_ARGS) $(RPI_ARGS) $(BMP_ARGS) -lm

all: codage libs

codage:	codage.c random.c turbocode.c rpi.c list.c ldpc.c basic.c
	gcc -o codage codage.c $(MAIN) $(MAIN_ARGS)

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