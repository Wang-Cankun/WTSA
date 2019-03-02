#ifndef _STRUCT_H
#define _STRUCT_H

#ifndef _GNU_SOURCE 
#define _GNU_SOURCE 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include <limits.h>
#include <sys/time.h>
#include <time.h>
#include "fib.h"

/***** Useful macros *****/

/* Compatibility of __attribute__ with non-GNU */
#ifndef __GNUC__
#  define __attribute__(x) /* Nothing */
#endif

/* Pretend that C has boolean type */
#define TRUE 1
#define FALSE 0
#define boolean unsigned char
#ifndef __cplusplus
#ifndef bool
#define bool unsigned char
#endif
#endif

#define MAX(a,b)  ((a)>(b)?(a):(b))
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define ABS(x)    ((x)>0?(x):-(x))

/* Variable and array allocation */
#define AllocVar(pt) (pt = xmalloc(sizeof(*pt)))
#define AllocArray(pt, size) (pt = xmalloc(sizeof(*pt) * (size)))
#define ReAllocArray(pt, size) (pt = xrealloc(pt, sizeof(*pt) * (size)))

/* Strings */
#define sameString(a, b) (strcmp((a), (b))==0)
/* Returns TRUE if two strings are same */

/* Debugging */
#define internalErr() errAbort("Internal error %s %d", __FILE__, __LINE__)

/* Constants */
#define LABEL_LEN 1000 
#define CLOSURE_SIZE 5000
#define MAX_SEQUENCE_LENGTH 5000
#ifndef NULL
#define NULL 0
#endif

/* Two major data types */
typedef float continuous;
typedef short discrete;

/* global data */
/*int extend_len;*/
continuous **arr;
discrete **arr_c;
discrete **arr_c1;
discrete **arr_v;
discrete **clo_matr;
discrete **clo_matr1;
char **sequences;
int  **height;
char *sequences_all;
char **sequences_regulon;
discrete size_closure;
continuous **p_markov;
continuous **td_markov;
/*simlation sequence - bingqiang*/
char **simulation_sequence;
discrete *symbols;
discrete* genes;
discrete* conds;
char **genes_n;
char **conds_n;
int sum_markov;
int all_id;
int* s_col;
int* promoter_length;
int oper_num_all;
int rows, cols, sigma,s_rows,s_cols,ver;
int clo_num;
char **locus_id;
bool *IsLengthEnough;
long sum_genome;
char **SequenceInfo;
int extend_len;
int is_extended; 

/*sequence weight*/
continuous *SequenceWeight;

/*we save the frequency matrix and score matrix so that we can reduce time complexity*/
discrete ***fre_all;
discrete ***fre2_all;
continuous ***sco_all;
continuous ***sco2_all;

/*improve the algorithm by change A T C G to 0 1 2 3*/
discrete **seq_matrix;
discrete ***fre_matrix;
discrete **height_matrix;


/*int clo_init_1, clo_init_2;*/
struct  dispos *profile1;

/***** Structures *****/

struct dyStack
/* dynamically allocated stack */
{
	int top;             /* top element index */
	int items[];		   /* data storage */
};
struct dyStack1
/* dynamically allocated stack */
{
	continuous top;             /* top element index */
	continuous items[];             /* data storage */
};

struct dispos
{
	continuous score;
	int x;
	int y;
	int z;
};

struct edgepos
{
	continuous score;
	int x;
	int y;
	int x1;
	int y1;
	int z;
};

/* edge between two genes */
typedef struct _Edge{
	int gene_one;
	int gene_two;
	int score;
}Edge;

/* biclustering block */
typedef struct _Block{
	struct dyStack *genes;
	struct dyStack *conds;
	int score;
	int oper_num;
	int block_rows;
	int block_cols;
	int block_rows_pre;
	double significance;
}Block;


typedef struct _Annotation{
        struct dyStack *init;
        struct dyStack *end;
        char **TF;
	int num;
}Annotation;

Annotation **anno;

typedef struct _Reference_genome{
        char **sequences_r;
	continuous **markov_matrix;
	char *oper_name;
	int seq_num;
	int *clo_num;
	continuous pvalue;
	int score;
	int motif_num;
	int motif_num_R;
}Reference_genome;

Reference_genome **genome;
int reference_oper_num;

typedef struct _Motif{
        struct dyStack *genes;
        struct dyStack *conds;
        float *score;
        int motif_number;
        double pvalue;
}Motif;


typedef struct _Closures{
        struct dyStack *sequence;
        struct dyStack *position;
	/*struct dyStack1 *score;*/
	continuous *score;
        long double  significance;
	long int pvalue;
	int closure_rows;
	int length;
	int size;
	char *name;
	continuous zscore;
        int zvalue;
        continuous enrich;
        int evalue;
	continuous motif_background_norm;
	continuous motif_known;
	discrete **frequency;
	continuous **scoreM;
	char **seed;
	/*struct dyStack *genes;
        struct dyStack *conds;
        int score;
        int block_rows;
        int block_cols;
        int block_rows_pre;
        double significance;*/
}Closures;

Closures **all;

/* holds running options */
typedef struct _Prog_options{
        char FN[LABEL_LEN];
	char BN[LABEL_LEN];
	char CN[LABEL_LEN];
	char GN[LABEL_LEN];
	char MN[LABEL_LEN];
	char HN[LABEL_LEN];
	char ZN[LABEL_LEN];
	bool IS_sequence;
	bool IS_SWITCH;
	bool IS_reference;
	bool IS_closure;
	bool IS_microarray;
	bool IS_reference_H;
	bool palindromic;
	bool mirror;
	bool IS_global;
	bool IS_local;
	bool zscore;
	bool expansion;
	bool ID;
	bool approximate;
	bool middle_enhance;
	bool no_enhance;
	bool FastVersion;
	bool SequenceWeight;
	bool checkPalin;
	bool RC;
	FILE* FP;
	FILE* FB;
	FILE* FC;
	FILE* FG;
	FILE* FM;
	FILE* FH;
	FILE* FZ;
	FILE* FJ;
	int COL_WIDTH;
	int SCH_BLOCK;
	int RPT_BLOCK;
	int MOTIFLENGTH;
	int FIRST;
	int SECOND;
	int THIRD;
	int TOPVERTICES;
	int threshold;
	int threshold_1;
	int threshold_2;
	int local2;
        int local3;
	int profile;
        int Low;
        int Up;
        int number;
	int simu;
	int range;
	int closure_enlarge_times;
	int conserve_threshold;
        int threshold_e2;
	int PalinLength;
	int Palin1;
	int Palin2;
	int Palin3;
        double end_weight;
	double FILTER;
	double TOLERANCE;
	double TOLERANCE1;
	double closure_threshold;
	continuous thre;
	continuous thre_pvalue;
	continuous zscore_thre;
	continuous conserve_background;
}Prog_options;

typedef unsigned short int bits16;
enum {UP=1, DOWN=2, IGNORE=3};

Prog_options* po;
/***** Helper functions *****/

void progress(char *format, ...)
/* Print progress message */
     __attribute__((format(printf, 1, 2)));

void verboseDot();
/* Print "i-am-alive" dot */

void err(char *format, ...)
/* Print error message but do not exit */
     __attribute__((format(printf, 1, 2)));

void errAbort(char *format, ...)
/* Print error message to stderr and exit */
     __attribute__((noreturn, format(printf, 1, 2)));

void uglyTime(char *label, ...);
/* Print label and how long it's been since last call.  Call with 
 * a NULL label to initialize. */

void *xmalloc ( int size );
/* Wrapper for memory allocations */

void *xrealloc ( void* ptr, int size );
/* Wrapper for memory re-allocations */

/* Stack-related operations */
struct dyStack *dsNew(int size);
struct dyStack1 *dsNew1(int size);
void dsPrint(struct dyStack *ds);
void dsPrint1(struct dyStack1 *ds);
void dsPush(struct dyStack *ds,int element);
void dsPush1(struct dyStack1 *ds,int element);
discrete** alloc2d(int rr, int cc);
continuous** alloc2dd(int rr, int cc);
char** alloc2c(int rr, int cc);
#define dsFree free
/* Release the stack data structure */

#define dsClear(pds) ((pds)->top = -1)
/* Remove all the data in the stack */

#define dsSize(pds) ((pds)->top + 1)
/* Return the size of the stack */

#define dsItem(pds, j) ((pds)->items[j])
/* Return the j-th item in the stack */

bool isInStack(struct dyStack *ds, int element);
int dsIntersect(struct dyStack *ds1, struct dyStack *ds2);

/* File-related operations */
void *addSuffix();
/* Return a string containing headsuffix */

FILE *mustOpen(const char *fileName, char *mode);
/* Open a file or die */

#endif
