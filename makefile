VER=2.0
DIST=WTSA$(VER)
PROGS=WTSA
SRCS=struct.c compare_sequence.c make_graph.c get_options.c fib.c write_file.c find_clique.c main.c closure_process.c pvalue.c read_file.c matrix_process.c 
OBJS=$(SRCS:.c=.o) 
CC=gcc


LDFLAGS=-static  -lm -fopenmp
CFLAGS=-O3 -Wall -ansi -I.  -DVER=$(VER) -fopenmp

all: $(PROGS)

WTSA: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)
.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -f $(PROGS)
	rm -f *.o
	rm -f *.closures

dist:
	$(MAKE) clean
	cd .. && tar czvf $(DIST).tar.gz $(DIST)/
test:
	$(MAKE)
	./BoBro -i example -l 14 -F
testP:
	$(MAKE)
	./BoBro -i example -P
testM:
	$(MAKE)
	./BoBro -i example -M
testG:
	$(MAKE)
	./BoBro -i example -G
testC:
	$(MAKE)
	./BoBro -i example -C
testa:
	$(MAKE)
	./BoBro -i example -j example_anno 
