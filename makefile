VER=3.0
DIST=BoBro$(VER)
PROGS=BoBro
SRCS=struct.c compare_sequence.c make_graph.c get_options.c fib.c write_file.c find_clique.c main.c closure_process.c pvalue.c read_file.c matrix_process.c 
OBJS=$(SRCS:.c=.o) 
CC=gcc -g


LDFLAGS= -static  -lm 
CFLAGS=-g -O0 -Wall -ansi -I.  -D VER=$(VER) 


all: $(PROGS)

BoBro: $(OBJS)
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
	./BoBro.exe -i test.fasta -l 14 -F
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
