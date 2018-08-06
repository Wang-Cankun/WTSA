#ifndef _CLOSURE_PROCESS_H
#define _CLOSURE_PROCESS_H

#include "struct.h"
continuous **d_markov (char **sequences_r, int seq_num);
continuous **markov (char **sequences_r, int seq_num);
continuous aver_score_closure(char **sequences_2, continuous **scoreM, continuous *score, int motif_number, int length_local_1);
int motif_num_closure (Reference_genome *genome_1, int motif_length, continuous **scoreM, continuous thre);
continuous *motif_num_R_closure (Reference_genome *genome_1, int motif_length, continuous **scoreM, continuous thre);
continuous get_pvalue(Reference_genome *genome_1, int motif_num, continuous *motif_num_R);
void print_operons (FILE *fw, char **sequences_regulon, Reference_genome **genome, int motif_num, int oper_num);
int block_cmpr_g(const void *a, const void *b);
void sort_genomes_list(Reference_genome **el, int n);
bool is_two_closures_regulon (Closures *a, Closures *b, Reference_genome *aa, Reference_genome *bb, int num_a, int num_b);


/*from closure_process*/
extern continuous **get_profile (char **sequence_temp, int a, int b, int motif_number);
extern continuous **impovre_profle (continuous **scoreM, double  length_local_1);
extern discrete **frequency_matrix (char **sequence_temp, int a, int b, int motif_number);
extern continuous **impovre_profle_palindromic (continuous **scoreM, double  length_local_1,  char **sequence_1, int a ,int motif_number);

#endif
