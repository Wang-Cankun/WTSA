#ifndef _CLUSTER_H
#define _CLUSTER_H

#include "struct.h"
/*prototypes*/
int cluster_1 (FILE *fw1,Edge **el, int n);
void print_bc_1 (FILE *fw1, Block* b, int num);
int post_processing_blocks_1 ( FILE *fw1, Block** bb, int num);

/*fron find_clique*/
extern void sort_block_list(Block **el, int n);
extern void block_init(Edge *e, Block *b,
                     struct dyStack *genes, struct dyStack *scores,
                     bool *candidates, const int cand_threshold,
                     int *components, struct dyStack *allincluster);
extern bits16 **profile;
extern void seed_current_modify (const discrete *s, bool *colcand, int* cnt, int components);
extern int intersect_row(const bool *colcand, discrete *g1, discrete *g2);
extern void scan_block (struct dyStack *gene_set, Block *b_ptr);
extern void seed_update (const discrete *s);

#endif

