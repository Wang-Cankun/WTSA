VER=2.0
DIST=WTSA$(VER)
PROGS=WTSA
SRCS=struct.c compare_sequence.c make_graph.c get_options.c fib.c write_file.c find_clique.c main.c closure_process.c pvalue.c read_file.c matrix_process.c 
OBJS=$(SRCS:.c=.o) 
CC=gcc -g


LDFLAGS= -static  -lm -fopenmp
CFLAGS=-g -O0 -Wall -ansi -I.  -DVER=$(VER) -fopenmp 


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
	./WTSA.exe -i test.fasta -l 14 -F
testP:
	$(MAKE)
	./WTSA -i example -P
testM:
	$(MAKE)
	./WTSA -i example -M
testG:
	$(MAKE)
	./WTSA -i example -G
testC:
	$(MAKE)
	./WTSA -i example -C
testa:
	$(MAKE)
	./WTSA -i example -j example_anno 
