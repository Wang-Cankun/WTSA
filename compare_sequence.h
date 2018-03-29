#ifndef _READ_ARRAY_H
#define _READ_ARRAY_H

/************************************************************************/
#include "struct.h"
/***********************************************************************/
/* prototypes  */
void init_dis();
void init_dis_1();
void parameter_choice();
static bool **matrix_no_continuous_equal(bool **matrix);
void compare_sequences(char **sequences);
static void pairwise_comparison_first ( bool **matrix3, int **matrix, bool *match, bool *match1, int lower, int upper, struct dispos *marray, struct edgepos *marray1, int mpos, int msize, int maxscore, int maxindex, int *startpos);
static void pairwise_comparison_second ( bool **matrix3, int **matrix, int **matrix1, bool *match, bool *match1, int lower, int upper, struct dispos *marray, struct edgepos *marray1, int mpos, int msize, int maxscore, int maxindex, int *startpos);
static void  pairwise_comparison_third (bool **matrix3, int **matrix, int **matrix1, bool *match, bool *match1, int lower, int upper, struct dispos *marray, struct edgepos *marray1, int mpos, int msize, int maxscore, int maxindex, int *startpos, int curelement1);
static void select_topvertex (int **matrix, bool **matrix2, bool **matrix3);
static void get_final_graph (int msize1, int **matrix,  bool **matrix2,struct edgepos *marray1);

/* from struct*/
extern continuous** alloc2dd(int rr, int cc);
extern discrete** alloc2d(int rr, int cc); 

/* from closure_process.c*/
continuous get_similarity_between_two_patterns(int seq1, int seq2, int pos1, int pos2, int motif_length);
continuous improve_similarity_between_two_patterns(int seq1, int seq2, int pos1, int pos2, int lower, int upper, continuous init);
/***********************************************************************/
#endif
