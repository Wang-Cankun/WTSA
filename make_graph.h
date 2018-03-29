#ifndef _MAKE_GRAPH_H
#define _MAKE_GRAPH_H

#include "struct.h"
/* prototypes */
int edge_cmpr(void *a, void *b);
void seed_update (const discrete *s);
int str_intersect_r (const discrete *s1, const discrete *s2);
void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min);
void fh_dump(struct fibheap *a, Edge **res);
void make_graph (const char *fn1);

/* global data */
Edge **edge_list;
Edge *edge_ptr;

/*from find_clique*/
extern int cluster (FILE *fw1,Edge **el, int n);

/*from matrix_process*/
extern int cluster_1 (FILE *fw1,Edge **el, int n);

/*from write_file*/
extern void print_params(FILE *fw1);
extern int report_closures(FILE *fw1, Closures **cc, int num_cc,Annotation **anno );
#endif

