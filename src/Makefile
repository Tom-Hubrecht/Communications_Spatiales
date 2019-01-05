LIB = "../lib/"
BIT_ARR_PATH = $(LIB)"BitArray"

all: codage packet

codage:	codage.c
	gcc -c codage.c -I $(BIT_ARR_PATH) -L $(BIT_ARR_PATH) -lbitarr

packet: packet.c
	gcc -c packet.c -I $(BIT_ARR_PATH) -L $(BIT_ARR_PATH) -lbitarr

libs:
	for libdir in $(LIB); do \
		$(MAKE) -C $$libdir; \
	done

clean:
	rm *.o

.PHONY: all clean libs