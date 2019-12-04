VER=2.0
DIST=Bobro$(VER)
PROGS=Bobro
SRCS=struct.c compare_sequence.c make_graph.c get_options.c fib.c write_file.c find_clique.c main.c closure_process.c pvalue.c read_file.c matrix_process.c 
OBJS=$(SRCS:.c=.o) 
CC=gcc -g


LDFLAGS= -static  -lm 
CFLAGS=-g -O0 -Wall -ansi -I.  -D VER=$(VER) 


all: $(PROGS)

Bobro: $(OBJS)
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
	./Bobro.exe -i test.fasta -l 14 -F
testP:
	$(MAKE)
	./Bobro -i example -P
testM:
	$(MAKE)
	./Bobro -i example -M
testG:
	$(MAKE)
	./Bobro -i example -G
testC:
	$(MAKE)
	./Bobro -i example -C
testa:
	$(MAKE)
	./Bobro -i example -j example_anno 
