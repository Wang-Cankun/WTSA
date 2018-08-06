#ifndef _CLOSURE_PROCESS_H
#define _CLOSURE_PROCESS_H

#include "struct.h"
/* prototypes  */
void read_annotation(FILE* fp);
void read_closures (FILE* fc);
void read_sequences(FILE* fp1 );
void read_reference_genome (FILE* fp);
void get_matrix_size (FILE* fp);
void read_labels (FILE* fp);
static int charset_add(discrete *ar, discrete s);
void init_dis_m();
void read_discrete (FILE* fp);
int scan_genome (continuous **scoreM, continuous AveScore, int motif_length, FILE* fp);
discrete *change_AGCT_to_num (char *seq, int length);

/*from closure_process*/
extern discrete **get_closure_matrix_1 (Closures **aa, int closure_id, continuous threshold);
extern continuous **get_palindromic_profile (continuous **scoreM, int length_local);
extern continuous palindromic_pro (char *consensus, int forepart);
extern bool palindromic_pro_1 (char *consensus, int forepart);

/*from pvalue*/
extern continuous **markov (char **sequences, int seq_num);
extern continuous **d_markov (char **sequences, int seq_num);

#endif
