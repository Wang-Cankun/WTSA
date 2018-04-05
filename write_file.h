#ifndef _WRITE_BLOCK_H
#define _WRITE_BLOCK_H

#include "struct.h"

/*prototypes*/
void print_params(FILE *fw1);
int report_closures(FILE *fw1, Closures **cc, int num_cc,Annotation **anno );
void sort_closures_list(Closures **el, int n);
int block_cmpr_1(const void *a, const void *b);
int report_regulon( FILE *fw,  Block** bb, int num);
void print_bc (FILE *fw1, Closures **cc, int num_cc, Block *b, int num);
void simu_markov(long double Motif_R_V[7], int seq_number,int length_local_1, continuous pp[5],continuous AveScore_V[7],continuous score_scan,continuous **scoreM);
static void print_regulon_horizonal( FILE *fw,  Block* bb, int num );
static void print_regulon_vertical( FILE *fw,  Block* bb, int num );


/* from find_clique */
extern void sort_block_list(Block **el, int n);
extern bits16 **profile;

/*from closure_process*/
extern continuous **get_profile (char **sequence_temp, int a, int b, int motif_number);
extern int get_num_TF (Closures *aa);
extern void print_frequency_anno (FILE *fw, char **aa, int num, int num1);
extern continuous **impovre_profle (continuous **scoreM, double length_local_1);
extern continuous **impovre_profle_palindromic (continuous **scoreM, double  length_local_1,  char **sequence_1, int a ,int motif_number);
extern continuous **impovre_profle_mirror (continuous **scoreM, double  length_local_1,  char **sequence_1, int a ,int motif_number);
extern int get_genome_num_from_closure (Closures *aa);
extern bool *clean_up_closures (Closures **aa, int closure_id, continuous threshold);
extern discrete **frequency_matrix (char **sequence_temp, int a, int b, int motif_number);
extern continuous **get_palindromic_profile (continuous **scoreM, int length_local);

/*from pvalue*/
extern continuous aver_score_closure(char **sequences_2, continuous **scoreM, continuous *score, int motif_number, int length_local_1);
extern void print_operons (FILE *fw, char **sequences_regulon, Reference_genome **genome, int motif_num, int oper_num);

/*from read_file*/
extern int scan_genome (continuous **scoreM, continuous AveScore, int motif_length, FILE* fp);
extern discrete *change_AGCT_to_num (char *seq, int length);
#endif
