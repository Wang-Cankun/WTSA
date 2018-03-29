#ifndef _CLUSTER_H
#define _CLUSTER_H

#include "struct.h"
/* prototypes */
int cluster (FILE *fw1,Edge **el, int n);
static int post_processing_blocks(FILE *fw1, Block **bb, int num);
void sort_block_list(Block **el, int n);
void block_init(Edge *e, Block *b,
                     struct dyStack *genes, struct dyStack *scores,
                     bool *candidates, const int cand_threshold,
                     int *components, struct dyStack *allincluster);
void seed_current_modify (const discrete *s, bool *colcand, int* cnt, int components);
static void update_colcand(bool *colcand,struct dyStack *cocluster,discrete *g1, discrete *g2);
int intersect_row(const bool *colcand, discrete *g1, discrete *g2);
int block_cmpr(const void *a, const void *b);
void scan_block (struct dyStack *gene_set, Block *b_ptr);
bits16 **profile;
void seed_update (const discrete *s);

/*from  struct */
extern void verboseDot();
discrete **alloc2d (int rr, int cc);
char **alloc2c (int rr, int cc);

/*from  write_file */
extern void print_params(FILE *fw1);
extern int report_closures(FILE *fw1, Closures **cc, int num_cc,Annotation **anno );
extern int report_regulon( FILE *fw,  Block** bb, int num);
extern void print_bc ( FILE *fw1, Closures **cc, int num_cc, Block *b, int num);
extern void sort_closures_list(Closures **el, int n);
extern int block_cmpr_1(const void *a, const void *b);
/*from closure_process*/
extern discrete **get_closure_matrix_1 (Closures **aa, int closure_id, continuous threshold);

/*from matrix_process*/
extern int post_processing_blocks_1 ( FILE *fw1, Block** bb, int num);
#endif

