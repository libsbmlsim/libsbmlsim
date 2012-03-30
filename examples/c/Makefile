PROG = test
OBJS = test.o
CC = gcc
CFLAGS = -Wall -I/usr/local/include -I../../src
LDFLAGS = -L/usr/local/lib -L../../build/src -lsbmlsim

.PHONY: all
all: $(PROG)

.SUFFIXES: .o .c
.c.o:
	$(CC) $(CFLAGS) -c $<

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean: 
	rm -f $(PROG) $(OBJS)
