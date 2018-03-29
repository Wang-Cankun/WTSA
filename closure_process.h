#ifndef _CLOSURE_PROCESS_H
#define _CLOSURE_PROCESS_H

#include "struct.h"
/*prototypes*/
discrete **frequency_matrix (char **sequence_temp, int a, int b, int motif_number);
continuous **get_profile (char **sequence_temp, int a, int b, int motif_number);
continuous **impovre_profle (continuous **scoreM, double length_local_1);
continuous get_IC (char **sequence_temp, int a, int b, int motif_number);
continuous IC_closures_1 (char **sequence_1, char **sequence_2, int a, int b, int motif_number_1, int motif_number_2, int p ,int q);
continuous IC_closures_2 (char **sequence_1, char **sequence_2, int a, int b, int motif_number_1, int motif_number_2,  int p ,int q);
char **get_sequences_from_closures (Closures *aa);
char **get_2L_sequeces_from_closures (Closures *aa);
continuous similarity_closures(Closures *aa, Closures *bb);
discrete **get_closure_matrix (Closures **aa, int closure_id, continuous threshold);
discrete **get_closure_matrix_1 (Closures **aa, int closure_id, continuous threshold);
void closure_clique (const char* fn);
int regulon (FILE *fw, Edge **el, int n);
int get_num_TF (Closures *aa);
void print_frequency_anno (FILE *fw, char **aa, int num, int num1);
continuous **impovre_profle_palindromic (continuous **scoreM, double  length_local_1,  char **sequence_1, int a ,int motif_number);
bool check_palindromic (char **sequence_palin, int a, int b, int motif_number,int forepart);
continuous palindromic_pro (char *consensus, int forepart);
bool palindromic_pro_1 (char *consensus, int forepart);
continuous mirror_pro (char *consensus, int forepart);
bool check_mirror (char **sequence_mirror, int a, int b, int motif_number,int forepart);
continuous **impovre_profle_mirror (continuous **scoreM, double  length_local_1,  char **sequence_1, int a ,int motif_number);
void update_matrix (continuous **matrix_score,discrete **matrix,int num);
int get_genome_num_from_closure (Closures *aa);
bool *clean_up_closures (Closures **aa, int closure_id, continuous threshold);
continuous *get_column_IC (Closures **aa, int clos);
continuous *get_column_IC_12 (Closures **aa, int clos_1, int clos_2, int startpos_1, int startpos_2);
continuous *get_column_IC_21 (Closures **aa, int clos_1, int clos_2, int startpos_1, int startpos_2);
continuous **get_palindromic_profle (continuous **scoreM, int length_local);
continuous get_similarity_between_two_patterns(int seq1, int seq2, int pos1, int pos2, int motif_length);
continuous improve_similarity_between_two_patterns(int seq1, int seq2, int pos1, int pos2, int lower, int upper, continuous init);
continuous **get_palindromic_profle (continuous **scoreM, int length_local);

/* global data */
int col_width;
Edge **edge_list;
Edge *edge_ptr;
extern bits16 **profile;

/* from make_graph*/
extern int edge_cmpr(void *a, void *b);
extern int str_intersect_r (const discrete *s1, const discrete *s2);
extern void seed_update (const discrete *s);
extern void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min);
extern void fh_dump(struct fibheap *a, Edge **res);

/* from find_clique*/
extern void block_init(Edge *e, Block *b, struct dyStack *genes, struct dyStack *scores,bool *candidates, const int cand_threshold, int *components, struct dyStack *allincluster);
extern void seed_current_modify (const discrete *s, bool *colcand, int* cnt, int components);
extern int intersect_row(const bool *colcand, discrete *g1, discrete *g2);
extern void scan_block (struct dyStack *gene_set, Block *b_ptr);

/*from  write_file */
extern int report_regulon( FILE *fw,  Block** bb, int num);

/*from pvalue*/
extern bool is_two_closures_regulon (Closures *a, Closures *b, Reference_genome *aa, Reference_genome *bb, int num_a, int num_b);
#endif
