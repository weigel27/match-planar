CC=g++
C=gcc

CFLAGS=-O2 -I. -march=native -msse3
CCFLAGS=-I. -O2 -march=native -msse3
LF = -lz
DEFS=

%.o: %.c
	$(C) -c $(CFLAGS) $(DEFS) $<

%.o: %.cc
	$(CC) -c $(CCFLAGS) $(DEFS) $<

TARGETS = match_ising_planar
OBJECTS = match_ising_planar.o blossom.o ising_options.o parse_opts.o ranvec.o zfstream.o

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJECTS)

match_ising_planar: match_ising_planar.o blossom.o ising_options.o parse_opts.o ranvec.o zfstream.o
	$(CC) -o match_ising_planar match_ising_planar.o blossom.o ising_options.o parse_opts.o ranvec.o zfstream.o $(LF)
