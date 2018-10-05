VER=1.0
DIST=wtsa$(VER)
PROGS=wtsa
SRCS=struct.c compare_sequence.c make_graph.c get_options.c fib.c write_file.c find_clique.c main.c closure_process.c pvalue.c read_file.c matrix_process.c 
OBJS=$(SRCS:.c=.o) 
CC=gcc -g


LDFLAGS= -static  -lm 
CFLAGS=-g -O0 -Wall -ansi -I.  -D VER=$(VER) 


all: $(PROGS)

wtsa: $(OBJS)
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
	./wtsa.exe -i test.fasta -l 14 -F
testP:
	$(MAKE)
	./wtsa -i example -P
testM:
	$(MAKE)
	./wtsa -i example -M
testG:
	$(MAKE)
	./wtsa -i example -G
testC:
	$(MAKE)
	./wtsa -i example -C
testa:
	$(MAKE)
	./wtsa -i example -j example_anno 
