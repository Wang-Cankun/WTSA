/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 16, 2010
 * Organize the parameters of wtsa
 */

/***************************************************************************/ 

#include "get_options.h"

/***************************************************************************/ 
/* ./wtsa or ./wtsa -h will pop up the usage of wtsa */
static const char USAGE[] = 
"===================================================================\n\
[Usage]\n\
$ ./wtsa -i filename [argument list]\n\
===================================================================\n\
[Input]\n\
-i : input file must be standard fasta file format\n\
[Optional Input]\n\
-j : input the TFBS annotation file\n\
-Z : input the background genome file\n\
===================================================================\n\
[Output]\n\
Commonly used parameters\n\
-l : motif length [5, ]\n\
     default: 14\n\
-L : minimum motif length\n\
     default: 14\n\
-U : maximum motif length\n\
     default: 14\n\
-o : number of closures to report (used under a specific input length)\n\
     default: 10\n\
-B : search reverse complement\n\
     default: TRUE\n\
Advanced parameters\n\
-n : top n closures under each length are used when L < U \n\
     default: 1\n\
-w : the weight of the two motif ends \n\
     default: 2\n\
-k : the minimum size of the initial motif seeds,\n\
     default: 3\n\
-c : consistency level of the motif seeds (0.5-1.0]\n\
     default: 1\n\
-s : the nunber of simulation times [5, ]\n\
     default: 5\n\
-u : the threshold of two closures' similarity socre (0,1]\n\
     default: 0.95\n\
-a : the upper limit of conservation level (N,10]\n\
     default: 9\n\
-N : the lower limit of conservation level (0,a)\n\
     default: 6\n\
-P : the flag of palindromic of TFBS\n\
     default: FALSE\n\
-M : the flag of mirror of TFBS\n\
     default: FALSE\n\
-G : the flag of global TF prediction\n\
     default: FALSE\n\
-C : the flag of local TF prediction\n\
     default: FALSE\n\
-E : the flag of expansion of closures base on the threshold 0.3-0.8\n\
     default: FALSE\n\
-A : the flag of approximation of pvalue calculation\n\
     default: FALSE\n\
-F : the flag of fast version of wtsa which just enhance two ends of motif\n\
     default: FALSE\n\
-W : the flag of considering sequences weight\n\
     default: FALSE\n\
-R : the range when we use [L,U]\n\
     default: 1\n\
-e : the times of seed alignments enlargement [1,3]\n\
     default: 3\n\
-b : the conserve level when search in background genome\n\
     default: 0.95\n\
===================================================================\n";

/***************************************************************************/ 
/* we assign initial values to the parameters */
static void init_options ()
{
	/* default parameters */
	strcpy(po->FN, " ");
	strcpy(po->BN, " ");
	strcpy(po->CN, " ");
	strcpy(po->GN, " ");
	strcpy(po->MN, " ");
	strcpy(po->HN, " ");
	strcpy(po->ZN, " ");
	po->COL_WIDTH = 3;
	po->TOLERANCE = 1;
    	po->FP = NULL;
	    po->FH =NULL;
	po->RPT_BLOCK = 10;
	po->SCH_BLOCK = 200;
        po->FIRST = 3; 
	po->SECOND = 4;
	po->THIRD = 3;
	po->FILTER =1;
	po->TOPVERTICES = 10;
	po->MOTIFLENGTH = 14;
	po->threshold = 5;
	po->threshold_1 = 5;
	po->threshold_2 = 5;
	po->local2 = 2; 
        po->local3 = 3;
        po->end_weight = 2;
	po->simu = 5;
        po->Low = 14;
        po->Up = 14;
	po->number = 1;
	po->IS_sequence = FALSE;
	po->IS_SWITCH = FALSE;
	po->IS_closure = FALSE;
	po->IS_reference = FALSE;
	po->IS_reference_H = FALSE;
	po->IS_microarray = FALSE;
	po->palindromic = FALSE;
	po->mirror = FALSE;
	po->IS_global = FALSE;
	po->IS_local = FALSE;
	po->expansion = FALSE;
	po->approximate = FALSE;
	po->ID =FALSE;
	po->middle_enhance = FALSE;
	po->FastVersion = FALSE;
	po->no_enhance = FALSE;
	po->SequenceWeight = FALSE;
	/*po->thre = 0.75;*/
	po->thre = 0.75;
	po->thre_pvalue = 0.05;
	/*po->closure_threshold = 0.7;*/
	po->closure_threshold = 0.7;
	po->range =1;
	po->closure_enlarge_times = 3;
	po->conserve_threshold = 9;
        po->threshold_e2 = 6;
	po->zscore = FALSE;
	po->conserve_background = 0.95;
	po->zscore_thre = 0;
	po->PalinLength = 30;
	po->Palin1 = 7;
	po->Palin2 = 5;
	po->Palin3 = 0;
	po->checkPalin = FALSE;
	po->RC = TRUE;
}
/***************************************************************************/ 
/* arrange input parameters */
void get_options (int argc, char* argv[])
{
	int op;
	bool is_valid = TRUE;
        AllocVar(po);
	/* get default value of parameters */
        init_options();
        while ((op = getopt(argc, argv, "i:j:r:g:m:FH:l:t:1:2:w:f:u:k:o:c:n:s:L:U:z:Z:d:b:p:e:a:N:PMGCEIA34R:Wq:v:x:y:Bh")) >0)
	{
		switch (op)
		{
                        case 'i': strcpy(po->FN, optarg); po->IS_sequence =TRUE; break;
			/*case '0': strcpy(po->FH, optarg); po->IS_height =TRUE; break;*/
			case 'j': strcpy(po->BN, optarg); po->IS_SWITCH = TRUE; break;
			case 'r': strcpy(po->CN, optarg); po->IS_closure = TRUE; break;
			case 'g': strcpy(po->GN, optarg);  po->IS_reference = TRUE;break;	  
			case 'm': strcpy(po->MN, optarg); po->IS_microarray = TRUE; break;
			case 'F': po->FastVersion = TRUE; break;
			case 'H': strcpy(po->HN, optarg); po->IS_reference_H = TRUE; break;
                        case 'l': po->MOTIFLENGTH = atoi(optarg); po->Low = atoi(optarg);po->Up = atoi(optarg); break;
			case 't': po->threshold = atoi(optarg); break;
			case '1': po->threshold_1 = atoi(optarg); break;
			case '2': po->threshold_2 = atoi(optarg); break;
                        case 'w': po->end_weight = atof(optarg); break;
			case 'f': po->FILTER = atof(optarg); break;
			case 'u': po->closure_threshold = atof(optarg); break;
			case 'k': po->COL_WIDTH = atoi(optarg); break;
			case 'c': po->TOLERANCE = atof(optarg); break;
			case 'n': po->number = atoi(optarg);break;
			case 'o':
				  po->RPT_BLOCK = atoi(optarg); 
			          po->SCH_BLOCK = 10 * po->RPT_BLOCK; break;
		        case 's': po->simu = atoi(optarg);break;
                        case 'L': po->Low = atoi(optarg); po->Up = po->Low; break;
                        case 'U': po->Up = atoi(optarg); break;
			case 'z': po->thre = atof(optarg); break;	  
			case 'Z': strcpy(po->ZN, optarg); po->zscore = TRUE; break;
			case 'd': po->zscore_thre = atof(optarg); break;	  
                        case 'b': po->conserve_background = atof(optarg); break;
			case 'p': po->thre_pvalue = atof(optarg); break;	  
			case 'e': po->closure_enlarge_times = atoi(optarg); break;	  
			case 'a': po->conserve_threshold = atoi(optarg); break;
                        case 'N': po->threshold_e2 = atoi(optarg); break;
			case 'P': po->palindromic = TRUE; break;	  
			case 'M': po->mirror = TRUE; break;	  
			case 'G': po->IS_global = TRUE; break;
			case 'C': po->IS_local = TRUE; break;
			case 'E': po->expansion = TRUE; break;
			case 'I': po->ID = TRUE; break;
			case 'A': po->approximate = TRUE; break;	  
			case '3': po->middle_enhance = TRUE; break;
			case '4': po->no_enhance = TRUE; break;
			case 'R': po->range = atoi(optarg); break;
			case 'W': po->SequenceWeight = TRUE; break;
			case 'q': po->PalinLength = atoi(optarg); po->checkPalin = TRUE; break;
			case 'v': po->Palin1 = atoi(optarg); break;
			case 'x': po->Palin2 = atoi(optarg); break;
			case 'y': po->Palin3 = atoi(optarg); break;
			case 'B': po->RC = FALSE; break;
			case 'h': puts(USAGE); exit(0); 
			default : is_valid = FALSE;
		}
	}
	/* ./wtsa will pop up the usage */
	if (argc==1) { puts(USAGE); exit(0);}
	/* read the input files  */
	if (po->IS_sequence) po->FP = mustOpen(po->FN, "r");
	/*if (po->IS_height) po->FH = mustOpen(po->FN, "r");*/
	if (po->IS_SWITCH)  po->FB = mustOpen(po->BN, "r");
	if (po->IS_closure)  po->FC = mustOpen(po->CN, "r");
	if (po->IS_reference) po->FG = mustOpen(po->GN, "r");
	if (po->IS_microarray) po->FM = mustOpen(po->MN, "r");
	if (po->IS_reference_H) po->FH = mustOpen(po->HN, "r");
	if (po->zscore) po->FZ = mustOpen(po->ZN, "r");
	/* option value range check */
	if ((po->TOLERANCE > 1) || (po->TOLERANCE <= .5))
	{
		err("-c conservativeness should be (.5,1]");
		is_valid = FALSE;
	}
	if (po->COL_WIDTH < 2 && po->COL_WIDTH != -1)
	{
		err("-k minimum column width should be >=2");
		is_valid = FALSE;
	}
        if ((po->MOTIFLENGTH < 5 && po->MOTIFLENGTH != -1 ))
        {
                err("-l minimum motif length should be >=5");
                is_valid = FALSE;
        }
        if (po->end_weight < 0)
        {
                err("-w the weight of the two ends should be >=0");
                is_valid = FALSE;
        }
        if (po->RPT_BLOCK <= 0)
	{
		err("-o number of blocks to report should be >0");
		is_valid = FALSE;
	}
        if (po->number <1)
	{
		err("-n number of closures of each lengeh should be >=1");
		is_valid = FALSE;
	}
	if (po->simu<5)
        {
                err("-s number of simulation should be >=5");
                is_valid = FALSE;
        }
        if (po->Low <= 4)
        {
                err("-L the lower length of motif should be >4");
                is_valid = FALSE;
        } 
        if (po->Up < po->Low)
        {
                err("-U the upper length of motif should be bigger than the lower length of motif");
                is_valid = FALSE;
        }
	if (po->closure_threshold <0 || po->closure_threshold >1)
	{
		err("-u the similarity cutoff between closures should be (0,1]");
		is_valid = FALSE;
	}
	if (po->closure_enlarge_times <= 0 || po->closure_enlarge_times >3)
	{
		err("-e the times of seed alignments enlargement should be [1,3]");
		is_valid = FALSE;
	}
	if (po->range <= 0)
	{
		err("-R the range when we use [L,U] should be [1,U-L]");
		is_valid = FALSE;
	}
	if (!is_valid) errAbort("Type -h to view possible options");
}

/***************************************************************************/ 
