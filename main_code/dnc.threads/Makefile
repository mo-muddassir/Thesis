PROGRAM = dnc
OBJECTS = dnc.o build_tree.o calc_forces.o particle_io.o move.o

CFLAGS = -Wall -oFast
LDLIBS = -lm -lpthread

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)

install:
	cp dnc ../bin/dnc

clean:
	rm -f $(PROGRAM) $(OBJECTS)
