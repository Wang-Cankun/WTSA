/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 12, 2010
 * Usage: This is part of BoBro package. Use, redistribution, modify without limitations

 * show how does the whole program work
 */

/***********************************************************************/

#include "main.h"

/***********************************************************************/

int main(int argc, char* argv[])
{
	/* Start the timer */
	uglyTime(NULL);
	/* get the program options defined in get_options.c */
	get_options(argc, argv);
	/* pop up the information of Bregulon */
	printf("\nWTSA %.2f: motif finding on genome scale (compiled "__DATE__" "__TIME__")\n\n", VER);
	/* read the fasta file and compare the sequences*/
	if (po->IS_sequence) {
		read_sequences(po -> FP);
		;
	}
	/* read the combined reference genome file if any (optional),
	 * so that we can refine and expand the predicted regulons by a pvalue threshold*/
	if (po->IS_reference) read_reference_genome (po->FG);
	if (po->IS_reference_H) read_reference_genome (po->FH);
	/* read the annotation file of input fasta file if any (optional),
	 * so that we can generate a STA line for each clousure */
	if (po->IS_SWITCH) {AllocArray (anno,s_rows);read_annotation(po -> FB);}
	/*intial the closures*/
 	int output = (po->Up - po->Low + 1)*3*po->RPT_BLOCK;
	AllocArray(all, output);
	/*identify the Motifs among different length*/
	int length_motif;
	if (po->FastVersion)
	{
		for (length_motif= po->Low; length_motif< (po->Up+1); length_motif += po->range, po->no_enhance = FALSE , po->middle_enhance = FALSE)
		{
			po->MOTIFLENGTH = length_motif;   
			compare_sequences(sequences);
	        init_dis();
			po->no_enhance = TRUE;
			extend_len=0;
			is_extended=0;
	        make_graph (addSuffix(po->FN, ".closures"));
			is_extended=1;
			make_graph (addSuffix(po->FN, ".closures"));
		}
	}
	else
	{
		for (length_motif= po->Low; length_motif< (po->Up+1); length_motif += po->range, po->no_enhance = FALSE , po->middle_enhance = FALSE)
		{
			extend_len=2;
			is_extended=1;
			po->MOTIFLENGTH = length_motif;   
			/*extend_len = 0;*/
			/* compare the input sequences in fasta format */
			compare_sequences(sequences);
			init_dis();
			/* find motif seeds from the graph we constructed*/
			make_graph (addSuffix(po->FN, ".closures"));
			
			po->middle_enhance = TRUE;
			compare_sequences(sequences);
			/*uglyTime("compare_sequences1", s_rows);*/
        	init_dis();
			/*uglyTime("init_dis1", s_rows);*/
            make_graph (addSuffix(po->FN, ".closures"));
			/*uglyTime("make_graph1", s_rows);*/

			po->no_enhance =TRUE;po->middle_enhance = FALSE;
			compare_sequences(sequences);
			/*uglyTime("compare_sequences2", s_rows);*/
            init_dis();
			/*uglyTime("init_dis2", s_rows);*/
	        make_graph (addSuffix(po->FN, ".closures"));
			/*extend_len = 1;
			init_dis();
			make_graph (addSuffix(po->FN, ".closures"));*/
			/*uglyTime("compare_sequences2", s_rows);*/
		}
	}
	/* find regulons base on the predicted closures */
	if (po->IS_reference_H)
	{
		init_dis_1();
		closure_clique (addSuffix(po->FN, ".regulons"));
	}
	/*if (po->IS_sequence) fclose(po->FP);*/
	if (po->IS_reference) fclose(po->FG);
	if (po->IS_SWITCH) fclose(po->FB);
	if (po->IS_reference_H) fclose(po->FH);
	free(po);
	return 0;
}
/***********************************************************************/
/*
for motif1, motif2 in all_motif_result
	for m1 in motif-1
		for m2 in motif-2
			if ( 1/2*motif_length < overlap < motif_length )
					m1_matrix = { min(start|end), max(start|end), m1_score + m2_score, length }
					or m1_matrix = { max(start|end), min(start|end), m1_score + m2_score, length } if ( m1_start > m2_start )
	new_motif_length = max_count_length in m1_matrix
	new_motif1 = {id, seq, m1_start, m1 + new_motif_length, new_motif_length_sequence}-**/