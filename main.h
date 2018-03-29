#ifndef _MAIN_H
#define _MAIN_H

#include "struct.h"
/* compare_sequence subroutine prototypes  */
extern void compare_sequences( char **sequences);
extern void read_sequences (FILE* fp1);
extern void read_annotation (FILE* fp);
extern void read_closures (FILE* fc);
extern void read_reference_genome (FILE* fp);
extern int get_options (int argc, char* argv[]);
extern void init_dis();
extern void init_dis_1();
extern void closure_clique (const char* fn);
extern void make_graph (const char* fn1);
extern void read_discrete (FILE* fp);
extern void get_matrix_size (FILE* fp);
extern void read_labels (FILE* fp);
#endif
