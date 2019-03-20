CC ?= gcc

PLATFORM := $(shell uname)
COMPILER := $(shell ($(CC) -v 2>&1) | tr A-Z a-z )

ifdef DEBUG
	OPT = -O0 -DDEBUG=1 --debug -g -ggdb
else
	ifneq (,$(findstring gcc,$(COMPILER)))
		OPT = -O4
		TGTFLAGS = -fwhole-program
	else
		OPT = -O3
	endif
endif

CFLAGS = -Wall -Wextra -I. $(OPT)
OBJFLAGS = -fPIC

all: libbmp.a

libbmp.o: libbmp.c libbmp.h

libbmp.a: libbmp.o
	ar -csru libbmp.a libbmp.o

%.o: %.c %.h
	$(CC) $(CFLAGS) $(OBJFLAGS) -c $< -o $@

clean:
	rm -rf libbmp.a *.o *.dSYM *.greg


# Comment this line out to keep .o files
.INTERMEDIATE: $(OBJS)
.PHONY: all clean test examples dev
