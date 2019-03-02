/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 15, 2010
 * the functions related to printing output to specific files
 */

#include "write_file.h"
#define INF 10000000
/*************************************************************************/
void print_params(FILE *fw1)
{
       /* fprintf(fw1, "#########################################\n#                                       #\n");
        fprintf(fw1, "#\tBoBro version %.2f output\t#\n", VER);
        fprintf(fw1, "#                                       #\n#########################################\n");
        fprintf(fw1, "\n****************************************\n");
        fprintf(fw1, "INPUT DATA SUMMARY\n");
        fprintf(fw1, "****************************************\n");
        fprintf(fw1,"\nDatafile: %s\n",po->FN);
        fprintf(fw1,"Alphabet: A G C T\n");
        fprintf(fw1,"Nucleotides composition:\tA: %.2f   G: %.2f   C: %.2f   T: %.2f\n",p_markov[1][0],p_markov[2][0],p_markov[3][0],p_markov[4][0]);
	fprintf(fw1,"Sequences number: %d\nNucleotides number: %d\n",s_rows,sum_markov);
        if (po->zscore)
	{ 
		fprintf(fw1,"Background Genome: %s\n",po->ZN);
		fprintf(fw1,"Background Genome Nuclectides number: %ld\n",sum_genome);
	}
        if (po->IS_SWITCH) fprintf (fw1,"Annotation data: %s\n",po->BN);
        if (po->IS_closure) fprintf (fw1,"Closures file: %s\n",po->CN);
	int i =0;
	fprintf(fw1,"\nSequence name\t\tLength\n");
	for (i=0;i<s_rows;i++)
	{
		/*fprintf (fw1,"%s\t\t%d\n",SequenceInfo[i],(s_col[i]+po->MOTIFLENGTH-1));*/
		/*fprintf (fw1,"%s\t\t%d\n",SequenceInfo[i],promoter_length[i]);
	}	
	
        fprintf(fw1, "\n****************************************\n");
        fprintf(fw1, "COMMAND LINE SUMMARY\n");
        fprintf(fw1, "****************************************\n");
        if (po->Low == po->Up)
	{
                fprintf(fw1, "\nParameters: -k %d; -c %.2f; -o  %d; -l %d; -u %2.2f; -e %d; -w %3.2f; -b %3.2f; -N %d;", po->COL_WIDTH, po->TOLERANCE, po->RPT_BLOCK,po->MOTIFLENGTH, po->closure_threshold, po->closure_enlarge_times, po->end_weight, po->conserve_background, po->threshold_e2);
		if (po->palindromic) fprintf(fw1, " -P");
		if (po->mirror) fprintf(fw1, " -M");
		if (po->IS_global) fprintf(fw1, " -G");
		if (po->IS_local) fprintf(fw1, " -C");
		if (po->IS_reference_H) fprintf(fw1, " -H");
		if (po->IS_microarray) fprintf(fw1, " -m");
		if (po->expansion) fprintf(fw1, " -E");
		if (SequenceWeight) fprintf(fw1, " -W");
		if (po->approximate) fprintf(fw1, " -A");
		if (po->FastVersion) fprintf(fw1, " -F");
                fprintf(fw1, "\n\nSeed size lower limit (-k)\t%d\nSeed conservation level (-c)\t%.2f\nOutput motif number (-o)\t%d\nMotif length (-l)\t%d\nMotif similarity lower limit (-u)\t%2.2f\nMotif refine times (-e)\t%d\nTwo ends weight (-w)\t%3.2f\nMotif conservation level in background (-b)\t%3.2f\nLower limit of motif conservation level (-N)\t%d\n", po->COL_WIDTH, po->TOLERANCE, po->RPT_BLOCK,po->MOTIFLENGTH, po->closure_threshold, po->closure_enlarge_times, po->end_weight, po->conserve_background, po->threshold_e2);
		if (po->palindromic) fprintf(fw1, "Panlindromic pattern of motif (-P)\tOn\n");
		if (po->mirror) fprintf(fw1, "Mirror pattern of motif (-M)\tOn\n");
		if (po->FastVersion) fprintf(fw1, "Fasta version (-F)\tOn\n");
	}
        else if (po->Low < po->Up)
	{
                fprintf(fw1, "\nParameters: -k %d; -c %.2f; -o %d; -n %d; -L %d; -U %d; -R %d; -u %2.2f; -e %d; -w %3.2f; -b %3.2f; -N %d;", po->COL_WIDTH, po->TOLERANCE, po->RPT_BLOCK, po->number, po->Low, po->Up, po->range, po->closure_threshold, po->closure_enlarge_times, po->end_weight, po->conserve_background, po->threshold_e2);
		if (po->palindromic) fprintf(fw1, " -P");
		if (po->mirror) fprintf(fw1, " -M");
                if (po->IS_global) fprintf(fw1, " -G");
                if (po->IS_local) fprintf(fw1, " -C");
                if (po->IS_reference_H) fprintf(fw1, " -H");
                if (po->IS_microarray) fprintf(fw1, " -m");
                if (po->expansion) fprintf(fw1, " -E");
		if (SequenceWeight) fprintf(fw1, " -W");
		if (po->approximate) fprintf(fw1, " -A");
		if (po->FastVersion) fprintf(fw1, " -F");
                fprintf(fw1, "\n\nSeed size lower limit (-k)\t%d\nSeed conservation level (-c)\t%.2f\nOutput motif number (-o)\t%d\nOutput motif number in each length (-n)\t%d\nMotif length lower limit (-L)\t%d\nMotif length upper limit (-U)\t%d\nMotif length step (-R)\t%d\nMotif similarity lower limit (-u)\t%2.2f\nMotif refine times (-e)\t%d\nTwo ends weight (-w)\t%3.2f\nMotif conservation level in background (-b)\t%3.2f\nLower limit of motif conservation level (-N)\t%d\n", po->COL_WIDTH, po->TOLERANCE, po->RPT_BLOCK, po->number, po->Low, po->Up, po->range, po->closure_threshold, po->closure_enlarge_times, po->end_weight, po->conserve_background, po->threshold_e2);
		if (po->palindromic) fprintf(fw1, "Panlindromic pattern of motif (-P)\tOn\n");
		if (po->mirror) fprintf(fw1, "Mirror pattern of motif (-M)\tOn\n");
		if (po->FastVersion) fprintf(fw1, "Fasta version (-F)\tOn\n");
	}
        fprintf(fw1, "\n");*/
}
/*************************************************************************/
int report_closures(FILE *fw1, Closures** cc, int num, Annotation** anno)
{
        sort_closures_list(cc, num);

	/*clean up the closures base on similarity scores*/
	bool *IS_duplicate;
	IS_duplicate = clean_up_closures(cc, num, po->closure_threshold);

        int closure_output=0, ii=0, jj=0,kk=0,extend_ii=0,extend_jj=0,extend_closure_output=0,total_instance=0;
        char **clo_TF;
	char **sequences_closures;
	/*num = MIN(num, po->RPT_BLOCK);*/
	/*fprintf (fw1,"Motif\tSeq\tstart\tend\tMotif\t\tScore\tInfo\n");*/
        while (ii< num && closure_output < po->RPT_BLOCK)
        {
                if (IS_duplicate[ii]) 
		{
			/*printf ("%d is duplicate\n", ii);*/
			ii++;
			continue;
		}
		if (cc[ii]->zscore < po->zscore_thre && po->zscore == TRUE)
		{
			ii++;
			continue;
		}
		closure_output++;
		/*fprintf (fw1,"\n\n*********************************************************\n");
                fprintf (fw1," Candidate Motif %3d\n",closure_output);
                fprintf (fw1,"*********************************************************\n\n");
                fprintf (fw1," Motif length: %d\n Motif number: %d\n Seed number: %d\n Motif Pvalue: %3.4LG (%ld)\n\n",cc[ii]->length+extend_len, cc[ii]->closure_rows, cc[ii]->size, cc[ii]->significance, cc[ii]->pvalue);
		*//*long double evalue;
		evalue = cc[ii]->significance * cc[ii]->closure_rows;
		fprintf (fw1," Motif Evalue: %3.4LG \n",evalue);*/
		
                /*if (po->zscore) fprintf (fw1," Seed Enrichment: %3.2f (%2.1f %2.5f)\n Seed Zscore: %3.2f\n",cc[ii]->enrich,cc[ii]->motif_known, cc[ii]->motif_background_norm, cc[ii]->zscore);
                fprintf (fw1,"\n------------------- Motif Seed------------------\n");
*/
                /*if (po->ID) fprintf (fw1,"#Seq\tposi\tID\tMotif\t\tScore\tInfo\tAnnotation\n");
		else fprintf (fw1,"#Seq\tposi\tMotif\t\tScore\tInfo\n");*/
		/*for (jj=0; jj<cc[ii]->size; jj++) 
		{
			for (kk=0;kk<cc[ii]->length+extend_len;kk++)
				fprintf (fw1, "%c", cc[ii]->seed[jj][kk]);
			fprintf (fw1,"\n");
		}	*/
                /*int i3 = get_num_TF (cc[ii]);
                clo_TF =alloc2c (i3,10);
                i3 = 0;
		sequences_closures = alloc2c(cc[ii]->closure_rows,cc[ii]->length);
                for (jj=0; jj<cc[ii]->size; jj++)
                {
                        int kkk=0;
			if (po->ID) fprintf (fw1,"%d\t%d\t%s\t",dsItem(cc[ii]->sequence,jj)+1,dsItem(cc[ii]->position,jj)+1,locus_id[dsItem(cc[ii]->sequence,jj)]);
			else fprintf (fw1,"%d\t%d\t",dsItem(cc[ii]->sequence,jj)+1,dsItem(cc[ii]->position,jj)+1);
                        for (kk=dsItem(cc[ii]->position,jj); kk< (dsItem(cc[ii]->position,jj)+cc[ii]->length); kk++)
			{
                                fprintf (fw1,"%c",sequences[dsItem(cc[ii]->sequence,jj)][kk]);
				sequences_closures[jj][kkk] = sequences[dsItem(cc[ii]->sequence,jj)][kk];	
				kkk++;
			}
                        fprintf (fw1,"\t%.2f",cc[ii]->score[jj]);
			fprintf (fw1,"\t%s",SequenceInfo[dsItem(cc[ii]->sequence,jj)]);
			
			if (po->IS_SWITCH)
                        {
                                int i1=0, i2=0;
                                for (kk = 0; kk<anno[dsItem(cc[ii]->sequence,jj)]->num;kk++)
                                {
                                        i1 = dsItem(anno[dsItem(cc[ii]->sequence,jj)]->init,kk);
                                        i2 = dsItem(anno[dsItem(cc[ii]->sequence,jj)]->end,kk);
                                        if ((dsItem(cc[ii]->position,jj) >= i1)&&(dsItem(cc[ii]->position,jj) <= i2))
                                        {
                                                clo_TF[i3++] = anno[dsItem(cc[ii]->sequence,jj)]->TF[kk];
                                                fprintf (fw1,"%s",anno[dsItem(cc[ii]->sequence,jj)]->TF[kk]);
                                                fprintf (fw1,"_%d",dsItem(anno[dsItem(cc[ii]->sequence,jj)]->init,kk));
                                                fprintf (fw1,"_%d ",dsItem(anno[dsItem(cc[ii]->sequence,jj)]->end,kk));
                                        }

                                }
                        }
                        fprintf (fw1,"\n");
                }*/
		
		/*fprintf (fw1,"\n------------------- Position weight matrix------------------\n");*/
		/*for (jj=1;jj<5;jj++)
		{
			if (jj==1) fprintf (fw1, "A\t");
			if (jj==2) fprintf (fw1, "G\t");
			if (jj==3) fprintf (fw1, "C\t");
			if (jj==4) fprintf (fw1, "T\t");
			for (kk=0;kk<cc[ii]->length;kk++)
			{
				fprintf (fw1, "%.2f\t",cc[ii]->scoreM[jj][kk]);
			}
			fprintf (fw1, "\n");
		}*/
		/*for (kk=0;kk<cc[ii]->length;kk++)
			fprintf(fw1,"\t%d",kk+1);
		fprintf (fw1,"\n");
		for (jj=1;jj<5;jj++)
		{
			if (jj==1) fprintf (fw1, "A\t");
			if (jj==2) fprintf (fw1, "G\t");
			if (jj==3) fprintf (fw1, "C\t");
			if (jj==4) fprintf (fw1, "T\t");
			for (kk=0;kk<cc[ii]->length;kk++)
			{
				fprintf (fw1, "%d\t",cc[ii]->frequency[jj][kk]);
			}
			fprintf (fw1, "\n");
		}
		fprintf (fw1,"IC");
		continuous temp_motif_num=0;
	       	continuous temp_ic=0;
		for (jj=1;jj<5;jj++)
			temp_motif_num+=cc[ii]->frequency[jj][1];
		for (kk=0;kk<cc[ii]->length;kk++)
                {
			temp_ic=0;
	                for (jj=1;jj<5;jj++)
			{
				temp_ic+=(cc[ii]->frequency[jj][kk]/temp_motif_num)*cc[ii]->scoreM[jj][kk];
			}
			fprintf (fw1, "\t%.2f",temp_ic);
                }
                fprintf (fw1, "\n");
		fprintf (fw1,"\n------------------- Consensus sequences------------------\n");
		char **consensus;
		int consensus_id=0, consensus_fre=0;
		consensus = alloc2c (5, cc[ii]->length+extend_len);
		for (kk=0;kk<cc[ii]->length+extend_len;kk++)
		{
			consensus_id=0;
			consensus_fre=0;
			for (jj=1;jj<5;jj++)
			{
				if (cc[ii]->frequency[jj][kk]>consensus_fre)
				{
					consensus_fre = cc[ii]->frequency[jj][kk];
					consensus_id = jj;	
				}
			}
			if (consensus_id==1) consensus[0][kk]='A';
			else if (consensus_id==2) consensus[0][kk]='G';
			else if (consensus_id==3) consensus[0][kk]='C';
			else if (consensus_id==4) consensus[0][kk]='T';
			if (consensus_fre>=0.75*cc[ii]->size) consensus[1][kk]=consensus[2][kk]=consensus[3][kk]=consensus[4][kk]='*';
			else if (consensus_fre>=0.5*cc[ii]->size) {consensus[1][kk]=consensus[2][kk]=consensus[3][kk]='*'; consensus[4][kk]=' ';}
			else if (consensus_fre>=0.25*cc[ii]->size) {consensus[1][kk]=consensus[2][kk]='*'; consensus[3][kk]=consensus[4][kk]=' ';}
			else if (consensus_fre>=0) {consensus[1][kk]='*'; consensus[2][kk]=consensus[3][kk]=consensus[4][kk]=' ';}
		}
		for (jj=0;jj<5;jj++)
		{
			for (kk=0;kk<cc[ii]->length+extend_len;kk++)
				fprintf (fw1, "%c", consensus[jj][kk]);
			fprintf (fw1,"\n");
		}	

                fprintf (fw1,"\n------------------- Aligned Motif ------------------\n");*/
                /*if (po->ID) fprintf (fw1,"Motif\tSeq\tposi\tID\tMotif\t\tScore\tInfo\tAnnotation\n");
		else fprintf (fw1,"Motif\tSeq\tstart\tend\tMotif\t\tScore\tInfo\n");*/

                int i3 = get_num_TF (cc[ii]);
                clo_TF =alloc2c (i3,10);
                i3 = 0;
		sequences_closures = alloc2c(cc[ii]->closure_rows,cc[ii]->length+extend_len);
                for (jj=0; jj<cc[ii]->closure_rows; jj++)
                {
                        int kkk=0;
			
			/*check whether positive or negative*/
			double positive=0, negative=0;
	                for (kk=dsItem(cc[ii]->position,jj); kk< (dsItem(cc[ii]->position,jj)+cc[ii]->length+extend_len); kk++)
        	        {
                                if (sequences[dsItem(cc[ii]->sequence,jj)][kk] == 'A' || sequences[dsItem(cc[ii]->sequence,jj)][kk] == 'a')
				{
					positive += cc[ii]->frequency[1][kkk]; 
					negative += cc[ii]->frequency[4][cc[ii]->length+extend_len-kkk-1];
				}
				else if (sequences[dsItem(cc[ii]->sequence,jj)][kk] == 'T' || sequences[dsItem(cc[ii]->sequence,jj)][kk] == 't')
				{
					positive += cc[ii]->frequency[4][kkk]; 
					negative += cc[ii]->frequency[1][cc[ii]->length+extend_len-kkk-1];
				}
				else if (sequences[dsItem(cc[ii]->sequence,jj)][kk] == 'G' || sequences[dsItem(cc[ii]->sequence,jj)][kk] == 'g')
				{
					positive += cc[ii]->frequency[2][kkk]; 
					negative += cc[ii]->frequency[3][cc[ii]->length+extend_len-kkk-1];
				}
				else
				{
					positive += cc[ii]->frequency[3][kkk]; 
					negative += cc[ii]->frequency[2][cc[ii]->length+extend_len-kkk-1];
				}
				kkk++;
        	        }
			
			if(is_extended == 0 && dsItem(cc[ii]->position,jj) > 0) {
				total_instance++;
			}
			kkk=0;
			if (positive >= negative)
			{
				if (po->ID) fprintf (fw1,">Motif-%d\t%d\t%d\t%d\t%s\t",closure_output,dsItem(cc[ii]->sequence,jj)+1,dsItem(cc[ii]->position,jj)+1,dsItem(cc[ii]->position,jj)+cc[ii]->length+extend_len, locus_id[dsItem(cc[ii]->sequence,jj)]);
				else fprintf (fw1,">Motif-%d\t%d\t%d\t%d\t",closure_output, dsItem(cc[ii]->sequence,jj)+1,dsItem(cc[ii]->position,jj)+1, dsItem(cc[ii]->position,jj)+cc[ii]->length+extend_len);
			}
			else
			{
				if (po->ID) fprintf (fw1,">Motif-%d\t%d\t%d\t%d\t%s\t",closure_output,dsItem(cc[ii]->sequence,jj)+1,dsItem(cc[ii]->position,jj)+cc[ii]->length+extend_len,dsItem(cc[ii]->position,jj)+1, locus_id[dsItem(cc[ii]->sequence,jj)]);
				else fprintf (fw1,">Motif-%d\t%d\t%d\t%d\t",closure_output, dsItem(cc[ii]->sequence,jj)+1,dsItem(cc[ii]->position,jj)+cc[ii]->length, dsItem(cc[ii]->position,jj)+1);
			}
                	for (kk=dsItem(cc[ii]->position,jj); kk< (dsItem(cc[ii]->position,jj)+cc[ii]->length+extend_len); kk++)
			{
				if (positive >= negative)
				{
                               		fprintf (fw1,"%c",sequences[dsItem(cc[ii]->sequence,jj)][kk]);
					sequences_closures[jj][kkk] = sequences[dsItem(cc[ii]->sequence,jj)][kk];	
				}
				else
				{
					int temp_posi = (dsItem(cc[ii]->position,jj)+cc[ii]->length+extend_len); 
                                	if (sequences[dsItem(cc[ii]->sequence,jj)][temp_posi-kkk-1] == 'A' || sequences[dsItem(cc[ii]->sequence,jj)][temp_posi-kkk-1] == 'a')
						fprintf (fw1,"T");
					else if (sequences[dsItem(cc[ii]->sequence,jj)][temp_posi-kkk-1] == 'T' || sequences[dsItem(cc[ii]->sequence,jj)][temp_posi-kkk-1] == 't')
						fprintf (fw1,"A");
					else if (sequences[dsItem(cc[ii]->sequence,jj)][temp_posi-kkk-1] == 'G' || sequences[dsItem(cc[ii]->sequence,jj)][temp_posi-kkk-1] == 'g')
						fprintf (fw1,"C");
					else
						fprintf (fw1,"G");
                	        	sequences_closures[jj][kkk] = sequences[dsItem(cc[ii]->sequence,jj)][temp_posi-kkk-1];
				
				}
				kkk++;
			}
                        fprintf (fw1,"\t%.2f",cc[ii]->score[jj]);
			fprintf (fw1,"\t%s",SequenceInfo[dsItem(cc[ii]->sequence,jj)]);
			
			if (po->IS_SWITCH)
                        {
                                int i1=0, i2=0;
                                for (kk = 0; kk<anno[dsItem(cc[ii]->sequence,jj)]->num;kk++)
                                {
                                        i1 = dsItem(anno[dsItem(cc[ii]->sequence,jj)]->init,kk);
                                        i2 = dsItem(anno[dsItem(cc[ii]->sequence,jj)]->end,kk);
                                        if ((dsItem(cc[ii]->position,jj) >= i1)&&(dsItem(cc[ii]->position,jj) <= i2))
                                        {
                                                clo_TF[i3++] = anno[dsItem(cc[ii]->sequence,jj)]->TF[kk];
                                                fprintf (fw1,"%s",anno[dsItem(cc[ii]->sequence,jj)]->TF[kk]);
                                                fprintf (fw1,"_%d",dsItem(anno[dsItem(cc[ii]->sequence,jj)]->init,kk));
                                                fprintf (fw1,"_%d ",dsItem(anno[dsItem(cc[ii]->sequence,jj)]->end,kk));
                                        }

                                }
                        }
                        fprintf (fw1,"\n");
                }
                if (po->IS_SWITCH)  print_frequency_anno (fw1, clo_TF, i3,cc[ii]->closure_rows );
		if (po->IS_reference_H)  
		{
			verboseDot();
			print_operons (fw1, sequences_closures, genome,cc[ii]->closure_rows,oper_num_all);
		}
		ii++;
                /*fprintf (fw1,"----------------------------------------------------\n");*/
        }
	int instance_length[total_instance];
	int instance_seq[total_instance];
	int instance_length_index = 0,compare_index=0,temp_overlap=0;
	if (is_extended == 0){
		while (extend_ii< num && extend_closure_output < po->RPT_BLOCK)
		{
			if (IS_duplicate[extend_ii]) 
			{
				printf ("%d is duplicate\n", ii);
				extend_ii++;
				continue;
			}
			if (cc[extend_ii]->zscore < po->zscore_thre && po->zscore == TRUE)
			{
				extend_ii++;
				continue;
			}
			for (extend_jj=0; extend_jj<cc[extend_ii]->closure_rows; extend_jj++)
			{
				if(dsItem(cc[extend_ii]->position,extend_jj) > 0)
				{
					instance_length[instance_length_index] = dsItem(cc[extend_ii]->position,extend_jj);
					instance_seq[instance_length_index] = dsItem(cc[ii]->sequence,jj);
					instance_length_index++;
					
				}
			}
			extend_closure_output++;
			extend_ii++;
		}
		for (instance_length_index=0; instance_length_index < total_instance-1;instance_length_index++)
		{
			for (compare_index=instance_length_index; compare_index < total_instance;compare_index++)
			{
				if(instance_seq[compare_index] == instance_seq[instance_length_index]) 
				{
					temp_overlap = abs(instance_length[compare_index] - instance_length[instance_length_index]);
					printf("%d\t\n",temp_overlap);
				}
				
			}
		}

	extend_len=2;


	}

        return closure_output;
}

/************************************************************************/
int block_cmpr_1(const void *a, const void *b)
/*compare function for qsort, decreasing by score*/
{
        if (po->zscore)
		return ((*(Closures **)b)->zvalue - (*(Closures **)a)->zvalue);
	else
		return ((*(Closures **)b)->pvalue - (*(Closures **)a)->pvalue);
}

void sort_closures_list(Closures **el, int n)
{
        qsort(el, n, sizeof *el, block_cmpr_1);
}
/************************************************************************/
int block_cmpr_oper(const void *a, const void *b)
/*compare function for qsort, decreasing by score*/
{
        return ((*(Block **)b)->oper_num - (*(Block **)a)->oper_num);
}

void sort_block_list_oper(Block **el, int n)
{
        qsort(el, n, sizeof *el, block_cmpr_oper);
}
/************************************************************************/
int report_regulon( FILE *fw,  Block** bb, int num)
{
	print_params(fw);
        fprintf (fw,"Input data: %s\n",po->FN);
        if (po->IS_SWITCH) fprintf (fw,"Annotation data: %s\n",po->BN);
        if (po->IS_closure) fprintf (fw,"Closures file: %s\n",po->CN);
        fprintf (fw,"Sequences number: %d\nNucleotides number: %d\n",s_rows,sum_markov);
        fprintf (fw,"Nucleotides composition:\tA: %.2f   G: %.2f   C: %.2f   T: %.2f\n",p_markov[1][0],p_markov[2][0],p_markov[3][0],p_markov[4][0]);
        /* if po->IS_reference is TRUE we should sort the regulons base on their containing operon numbers */
	if (po->IS_reference) sort_block_list_oper (bb,num);
	else sort_block_list(bb, num);

        int i, j,k;
        int n = MIN(num, po->RPT_BLOCK);
        bool flag;
        Block **output;
        AllocArray(output, n);
        Block **bb_ptr = output;
        Closures **cc;
        AllocArray(cc, n);
        Block *b_ptr;
        int cur_rows, cur_cols;
        int inter_rows, inter_cols;
        /* the major post-processing here, filter overlapping blocks*/
        i = 0; j = 0;
        while (i < num && j < n)
        {
                b_ptr = bb[i];
                cur_rows = b_ptr->block_rows;
                cur_cols = b_ptr->block_cols;
                flag = TRUE;
                k = 0;
                while (k < j)
                {
                        inter_rows = dsIntersect(output[k]->genes, b_ptr->genes);
                        inter_cols = dsIntersect(output[k]->conds, b_ptr->conds);
                        if (inter_rows*inter_cols > po->FILTER*cur_rows*cur_cols) flag = FALSE; break;
                        k++;
                }
                i++;
		if (flag)
                {
/*        printf ("%d\t%d\n",dsItem(all[dsItem(b_ptr->genes,0)]->sequence,0),dsItem(all[dsItem(b_ptr->genes,0)]->position,0));*/
                        if (po->IS_reference) print_regulon_vertical (fw,b_ptr,j++);
			else	print_regulon_horizonal (fw,b_ptr,j++);
                        *bb_ptr++ = b_ptr;
			/*verboseDot();*/
                }
        }
        return j;
}

/************************************************************************/
/* Identified clusters are backtraced to the original data, by putting the clustered vectors together, identify common column */
void print_bc (FILE *fw1, Closures **cc, int num_cc, Block* b, int num)
{
        int i;
        int block_rows;
	int background_num;
	continuous motif_background=0, motif_known=0, enrichment=0, zscore=0, motif_background_norm=0;

        block_rows = b->block_rows;
        char **sequences_1;
        sequences_1=(char**)malloc(sizeof(char*)*block_rows);
        for(i=0;i<b->score;i++)
                sequences_1[i]=(char*)malloc(sizeof(char)*50);
        char *seq;
        AllocArray(seq,100);
        int *pos,*dd,*dd1;
        continuous  *IC, *IC_sum, *IC_sum1;
	continuous AveScore=0;
        AllocArray (pos,100);
        AllocArray (dd,100);
        AllocArray (IC,4);
        AllocArray (IC_sum,100);
        AllocArray (IC_sum1,20);
        AllocArray (dd1,100);
        int j,k=0,d=0,d1,left,right;
        for(i=0;i<100;i++)
                dd[i]=0;
	/*print the profile corresponding to the accompany matrix*/
        for (i=0;i<b->score;i++)
        {
                d = profile1[genes[dsItem(b->genes, i)]].x;
                d1 = profile1[genes[dsItem(b->genes, i)]].y;
                k=0;
                if (i==0) dd1[0]=d1;
                else
                {
                        if (dsItem(b->genes, 0) <= dsItem(b->genes, i))
                                dd1[i]=d1+arr_c1[dsItem(b->genes, 0)][dsItem(b->genes, i)];
                        else
                                dd1[i]=d1-arr_c1[dsItem(b->genes, 0)][dsItem(b->genes, i)];
                }
                left = MAX(0,dd1[i]-4);
                int length1=strlen(sequences[d])-1;
                right = MIN(length1,dd1[i]+po->MOTIFLENGTH+4);
                if (dd1[i]<4)
                        for (j=0;j<4-dd1[i];j++)
                        {
                                sequences_1[i][k]='N';
                                k++;
                        }
                for (j=left;j<right;j++)
                {
                        sequences_1[i][k]=sequences[d][j];
                        k++;
                }
                if (length1<dd1[i]+po->MOTIFLENGTH+4)
                        for (j=0;j<(dd1[i]+po->MOTIFLENGTH+4-length1);j++)
                        {
                                sequences_1[i][k]='N';
                                k++;
                        }
                dd1[i]=dd1[i]-4;
        }
	int num1=0,num11=0;
        double num4=0, average1=0,average=0;
        IC[0]=IC[1]=IC[2]=IC[3]=0;
        for (i=0;i< (po->MOTIFLENGTH+8);i++)
        {
                IC[0]=IC[1]=IC[2]=IC[3]=0;
                IC_sum[i]=0;
                for (j=0;j<block_rows;j++)
                {
                        if (sequences_1[j][i]=='A') IC[0]++;
                        if (sequences_1[j][i]=='T') IC[1]++;
                        if (sequences_1[j][i]=='C') IC[2]++;
                        if (sequences_1[j][i]=='G') IC[3]++;
                }
                for (k=0;k<4;k++)
                        if(IC[k]>0)
                        {
                                IC[k]=(IC[k]/block_rows)*log(4*IC[k]/block_rows);
                                IC_sum[i]+=IC[k];
                        }
        }
	/*local optimal of IC*/
        double length, length_local=0,two_end=0;
        left= MAX (po->MOTIFLENGTH,6);
        for (length=left;length<po->MOTIFLENGTH+1;length++)
        {
                for (i=0;i<po->MOTIFLENGTH+8-length;i++)
                {
                        IC_sum1[i]=0;
                        two_end=0;
                        IC_sum1[i]+=po->end_weight*IC_sum[i];
                        two_end+=IC_sum[i];
                        for (j=i+1;j<i+floor(length/3);j++)
                        {
                                IC_sum1[i]+=po->end_weight*IC_sum[j];
                                two_end+=IC_sum[j];
                        }
                        for (j=i+floor(length/3);j<i+length-ceil(length/3);j++)
                                IC_sum1[i]+=IC_sum[j];
                        for (j=i+length-ceil(length/3);j<i+length;j++)
                        {
                                IC_sum1[i]+=po->end_weight*IC_sum[j];
                                two_end+=IC_sum[j];
                        }
                        if (IC_sum1[i] >= num4)
                        {
                                num4=IC_sum1[i];
                                average1=two_end/(floor(length/3)+ceil(length/3));
                                num11=i;
                        }
                }
                if (average1>=average)
                {
                        average=average1;
                        num1=num11;
                        length_local=length;
                }
        }
	/*the information which bingqiang will need*/
        for (i=0; i<dsSize(b->genes); i++)
                d = profile1[genes[dsItem(b->genes, i)]].x;

	/***********************get closures for seed**********************/
        Closures *cctemp;
        AllocVar(cctemp);
        int ii,t;        
        int  seq_number = s_rows; 
        int length_local_1 = length_local+extend_len, Motif_Scan_V[7];
        /*continuous Motif_R_V[7],pvalue_V[7];*/
	long double Motif_R_V[7],pvalue_V[7];
        continuous pp[5],AveScore_V[7],score_scan=0, score_scan_RC=0;
        discrete **frequency;
        continuous **scoreM;
        continuous **scoreM_RC;
        Motif* c;
        AllocVar(c);
        struct dyStack *gene, *cond;
        struct dyStack *gene2, *cond2 ;
        continuous score[CLOSURE_SIZE];
        int motif_number= b->score,motif_number_2=0;
        long double  pvalue=1; pvalue=1;
        gene = dsNew(CLOSURE_SIZE);
        cond = dsNew(CLOSURE_SIZE);
        gene2 = dsNew(CLOSURE_SIZE);
        cond2 = dsNew(CLOSURE_SIZE);
	if (num==0)  size_closure = motif_number;
        for (i=0; i<motif_number; i++)
	{
                dsPush(gene,dsItem(b->genes,i));
                dsPush(cond,dsItem(b->conds,i));
                score[i]=0;
        }

	cctemp->frequency= alloc2d (5,length_local_1);
	cctemp->scoreM = alloc2dd (5,length_local_1);
        frequency = alloc2d (5,length_local_1); 
	scoreM = alloc2dd (5,length_local_1);
	scoreM_RC = alloc2dd (5,length_local_1);
	cctemp->seed = alloc2c (motif_number, length_local_1);

        for (i=0; i<5; i++)
	{
		for (j=0; j<length_local_1; j++)
		{
			scoreM[i][j]=0;
			scoreM_RC[i][j]=0;
		}
	}

        /* 1000: 3 times closure improvement begion : ii */
        for (ii=1;ii< po->closure_enlarge_times+1;ii++)
	{
                /*1100: load current seed or motifs in sequences_2: begin*/
                char **sequences_2;
                sequences_2 = alloc2c (motif_number, length_local_1);
                if (ii==1)
                        for(i=0;i<motif_number;i++)
                                for (j=0; j<length_local_1; j++)
                                        sequences_2[i][j] = sequences_1[i][num1+j];
                else
		{
                        motif_number=motif_number_2;
                        for(i=0;i<motif_number;i++)
                                for (j=0; j<length_local_1; j++)
                                        sequences_2[i][j]=sequences[dsItem(gene2,i)][dsItem(cond2,i)+j];
                }         
                /*1200: calculate the frequency matrix and log matrix begin*/
                scoreM = get_profile (sequences_2, 5 , length_local_1, motif_number);
		if (ii==1) 
		{
			frequency = frequency_matrix (sequences_2, 5, length_local_1, motif_number);
			cctemp->frequency = frequency;
			for (j=1;j<5;j++)
			{
				for (k=0;k<length_local_1;k++)
				{
					cctemp->scoreM[j][k] = scoreM[j][k];
				}
			}
			cctemp->seed = sequences_2;
		}
                /* 1300: improve the scoring log matrix begin*/
		if (po->palindromic)	scoreM = impovre_profle_palindromic (scoreM, length_local_1,  sequences_2, 5 ,motif_number);
		else if (po->mirror)    scoreM = impovre_profle_mirror (scoreM, length_local_1,  sequences_2, 5 ,motif_number);
		else 
		{	
			scoreM = impovre_profle (scoreM, length_local_1);
                	scoreM_RC = get_palindromic_profile (scoreM , length_local_1);
		}

                /*1400: scoring motifs and calculate threshod for closure begin*/
		AveScore = aver_score_closure(sequences_2,  scoreM, score, motif_number, length_local_1);
/*		printf ("%3.2f\n",AveScore);*/
                for (i=0;i<8;i++)
		{
                        AveScore_V[i]=(0.3+0.1*i)*AveScore;
                        Motif_Scan_V[i]=0; 
			Motif_R_V[i]=0;
                }
                /*1500: scan data to built closures begin*/
                gene = dsNew(CLOSURE_SIZE);
                cond = dsNew(CLOSURE_SIZE);
                motif_number=0;
		for (i=0;i<seq_number;i++)
                {
			if (po->RC)
			{
			        for (j=0;j<strlen(sequences[i])-length_local_1+1;j++)   
				{
                	                score_scan=0;
                        	        score_scan_RC=0;
                                	for (k=0;k<length_local_1;k++)
					{
						score_scan += scoreM[seq_matrix[i][k+j]+1][k];
						score_scan_RC += scoreM_RC[seq_matrix[i][k+j]+1][k];
					}
					if ((score_scan>=AveScore_V[po->threshold_e2-3] || score_scan_RC >= AveScore_V[po->threshold_e2-3])&& motif_number < CLOSURE_SIZE-10)
					{
        	                                for (k=6;k>=0;k--)
                	                                if (score_scan>AveScore_V[k] || score_scan_RC > AveScore_V[k])
                        	                                Motif_Scan_V[k]++;
                                	        dsPush(gene,i);
                                        	dsPush(cond,j);
	                                        score[motif_number]=MAX(score_scan,score_scan_RC);
        	                                motif_number++;
                	                }
                        	}
			}
			else
			{
			        for (j=0;j<strlen(sequences[i])-length_local_1+1;j++)   
				{
                	                score_scan=0;
                                	for (k=0;k<length_local_1;k++)
						score_scan += scoreM[seq_matrix[i][k+j]+1][k];
					if (score_scan>=AveScore_V[po->threshold_e2-3] && motif_number < CLOSURE_SIZE-10)
					{
        	                                for (k=6;k>=0;k--)
                	                                if (score_scan>AveScore_V[k])
                        	                                Motif_Scan_V[k]++;
                                	        dsPush(gene,i);
                                        	dsPush(cond,j);
	                                        score[motif_number]=score_scan;
        	                                motif_number++;
                	                }
                        	}
			}
				
		}
                /*1600: for first two steps of closure-build only begin*/
		if (ii < po->closure_enlarge_times)
		{
                        gene2 = dsNew(CLOSURE_SIZE);
                        cond2 = dsNew(CLOSURE_SIZE);
                        if (motif_number<size_closure)
			{
                                motif_number_2=motif_number;
                                for (i=0;i<motif_number;i++)
				{
                                        dsPush(gene2,dsItem(gene,i));
                                        dsPush(cond2,dsItem(cond,i));
                                }
                        }
			else
			{
                                int i_max_1=0; continuous score_max_1=0;
                                motif_number_2=size_closure;
                                for (i=0;i<motif_number_2;i++)
				{
                                        i_max_1=0;score_max_1=0;
                                        for (j=0; j<motif_number; j++)
                                                if (score[j]>score_max_1) 
						{
                                                        i_max_1=j;
                                                        score_max_1=score[j];
                                                }
                                        dsPush(gene2,dsItem(gene,i_max_1));
                                        dsPush(cond2,dsItem(cond,i_max_1));
                                        score[i_max_1]=0;
                                }
                        }
                }
        }

		/*scan the background and get the motif number in whole gemone*/
        if (po->zscore)
        {
                /*continuous zscore,enrichment,motif_background,motif_known;*/
                motif_known = scan_genome (scoreM, AveScore, length_local_1, po->FP);
                background_num = scan_genome (scoreM, AveScore, length_local_1, po->FZ);
                if (background_num==0)
                {
                        printf("Do not find any instances in background\n");
                        enrichment = INF; zscore = INF;
                }
                else
                {
                        motif_background = background_num;
			motif_background_norm = (sum_markov*motif_background)/sum_genome;
                        /*printf ("%ld\t%d\t%.2f\t%.2f\t%.2f\n",sum_genome,sum_markov,motif_background,motif_known,motif_background_norm);*/
                        enrichment = (motif_known*sum_genome)/(sum_markov*motif_background);
                        zscore = (motif_known-sum_markov*motif_background/sum_genome)/sqrt(sum_markov*motif_background/sum_genome);                               
                }
                /*printf ("%d\t%d\t%ld\t%d\t%.2f\t%.2f\n",background_num,motif_number,sum_genome,sum_markov,enrichment,zscore);*/
        }
        /* 1000: 3 times closure improvement end : ii */

        /* 2000:  simulation on markov data begin */
	int length_ave=0,numberfore=1,number_pre=1,number_now=1,randomNum,simulation=po->simu,Rt,num_all_R=0;
        continuous length_ave_1=0;
        for (i=0;i<seq_number;i++)
                length_ave_1=length_ave_1+strlen(sequences[i]);
        length_ave=ceil(length_ave_1/seq_number); 

        srand((unsigned)time(NULL));
        char *randomdata;
        AllocArray (randomdata,length_ave);
	discrete *randomdata_number;
        randomdata_number = change_AGCT_to_num(randomdata,length_ave);
        for (Rt=0;Rt<4000*simulation;Rt++)
	{
                for (i=0;i<4;i++)
		{
                        for(j=0;j<length_ave;j++)  
			{
                                num_all_R++;
                                for (k=1;k<5;k++)
				{
                                        pp[k]=td_markov[numberfore][k];
					number_pre = number_now;
				}
                                randomNum=rand()%100;
				/*for 3-order markov, pre = (pre - 1) * 4 + now */
                                if (randomNum<pp[1]*100) 
				{
					randomdata[j]='A';numberfore=(number_pre-1)*4+1;
					number_now = 1;
				}
                                else if (randomNum<(pp[1]+pp[2])*100) 
				{
					randomdata[j]='G';numberfore=(number_pre-1)*4+2;
					number_now = 2;
				}
                                else if (randomNum<(pp[1]+pp[2]+pp[3])*100) 
				{
					randomdata[j]='C';numberfore=(number_pre-1)*4+3;
					number_now = 3;
				}
                                else 
				{
					randomdata[j]='T';numberfore=(number_pre-1)*4+4;
					number_now = 4;
				}
                        } 
			randomdata_number = change_AGCT_to_num(randomdata,length_ave);
                        for(j=0;j<length_ave-length_local_1+1;j++)  
			{
                                score_scan=0;
				for (k=0;k<length_local_1;k++)
					score_scan += scoreM[randomdata_number[k+j]+1][k];
				for (k=6;k>=0;k--)
                                        if (score_scan>AveScore_V[k])
                                                Motif_R_V[k]=Motif_R_V[k]+1;
                        } 
                }
	}
	for (t=6;t>=0;t--) 
                Motif_R_V[t]=(Motif_R_V[t]*seq_number)/(4*(Rt));

        /* 2000:  simulation on markov data end   */

        /* 3000:  pvalue calculating begin   */     /*note: 4 corresponding to 0.7; 3 corresponding to 0.6; ...*/
        /*continuous Scan_tmp,R_tmp,poisson,one_tmp=1,iiii;*/
        long double Scan_tmp,R_tmp,poisson,one_tmp=1,iiii;
	int tt=0;
        for (t=po->conserve_threshold-3;t>=po->threshold_e2-3;t--)
	{
                pvalue_V[t]=0;
                Scan_tmp=Motif_Scan_V[t];
                R_tmp=Motif_R_V[t];
                while (R_tmp >600) 
		{        
			Scan_tmp=Scan_tmp/2;
			R_tmp=R_tmp/2;
		}
		if (R_tmp == 0) R_tmp++;
/*		printf ("%3.2f\n",R_tmp);*/
                poisson=one_tmp/exp(R_tmp);
		if (po->approximate)
			j=1;
		else 
			j=300;
                for (iiii=0;iiii<(int)(Scan_tmp)+j;iiii++)
		{
                        if (iiii>(int)(Scan_tmp)-1) pvalue_V[t]=pvalue_V[t]+poisson;
                        poisson=poisson*R_tmp/(iiii+1);
                }
                if (pvalue>pvalue_V[t])
		{
                        pvalue=pvalue_V[t];
			tt = t;
		}
        }
/*	printf ("%LG\t%LG\t%LG\t%LG\n",pvalue,pvalue_V[tt],poisson,R_tmp);*/
        /*adjust the output closures because the closure corresponding to the minimal pvalue will be too small*/
	if (po->expansion)
	{
		tt=0;
	        for (t=0;t<6;t++)
                	if ((Motif_Scan_V[t]-Motif_Scan_V[t+1])/(Motif_R_V[t]-Motif_R_V[t+1]) > 1.2)
				{tt=t;	break;}
	}
        /* 3000:  pvalue calculating end   */
	/* 4000:  fix the final closure and write in struct Closure begin  */
        struct dyStack *gene1, *cond1 ;
        continuous score1[CLOSURE_SIZE];
        gene1 = dsNew(CLOSURE_SIZE);
        cond1 = dsNew(CLOSURE_SIZE);
        int i_max=0;
        continuous score_max=AveScore_V[tt];
        int motif_number_1=0;
        for (i=0;i<Motif_Scan_V[tt];i++)
	{
                i_max=0;score_max=AveScore_V[tt];
                for (j=0;j<motif_number;j++)
                        if (score[j]>score_max)
			{
                                score_max=score[j];
                                i_max=j;
                        }
                dsPush(gene1,dsItem(gene,i_max));
                dsPush(cond1,dsItem(cond,i_max));
                score1[motif_number_1]=score[i_max];
                motif_number_1++;
                score[i_max]=0;
        }/*sort and cut closure end*/

        /*Closures *cctemp;
        AllocVar(cctemp);*/
        cctemp->sequence = dsNew(motif_number_1);
        cctemp->position = dsNew(motif_number_1);
        cctemp->score = (continuous*)malloc(sizeof(continuous)*motif_number_1);
	cctemp->motif_background_norm = motif_background_norm;
	cctemp->motif_known = motif_known;
        int kk;
        for ( kk=0; kk < motif_number_1; kk++)
	{
                dsPush(cctemp->sequence,dsItem(gene1,kk));
                dsPush(cctemp->position,dsItem(cond1,kk));
                cctemp->score[kk] = score1[kk];
                cctemp->closure_rows = motif_number_1;
                cctemp->significance = pvalue;
                cctemp->pvalue = -(100*log10(pvalue));
                cctemp->length = po->MOTIFLENGTH;
		cctemp->name = po->FN;
		cctemp->size = b->score;
		cctemp->zscore = zscore;
                cctemp->zvalue = 100*zscore;
                cctemp->enrich = enrichment;
                cctemp->evalue = 100*enrichment;

        }
        cc[num_cc] = cctemp;
}
/******************************************************************/
static void print_regulon_horizonal( FILE *fw,  Block* bb, int num )
{
	/*printf ("%d\t%d\n",dsItem(all[dsItem(bb->genes,0)]->sequence,0),dsItem(all[dsItem(bb->genes,0)]->position,0));*/
	int i,j,kk=0,kk_1;
	long int k=0;
	verboseDot();
	fprintf (fw,"\n\n*********************************************************\n");
	fprintf (fw,"Candidate Regulon %d:\t ",(num+1));
	for (i=0;i< (bb->block_rows); i++)
		 fprintf (fw,"closure %d ",(dsItem(bb->genes,i)+1));
	fprintf (fw, "\n");
	fprintf (fw,"*********************************************************\n\n");
	
	for (i=0;i<(bb->block_rows);i++)
		kk += all[dsItem(bb->genes,i)]->closure_rows;
	/*delete the redundancy*********************************************/	
	bool *isclo;
	AllocArray(isclo, s_rows*strlen(sequences[0]));
	for (i=0;i<s_rows*strlen(sequences[0]);i++)
		isclo[i] = FALSE;	
	bool *isregu;
	AllocArray(isregu,s_rows);
	for (i=0;i<s_rows;i++)
		isregu[i] = FALSE;
	kk = 0,kk_1=0;
	for (i=0;i<(bb->block_rows);i++)	
	{
		for (j=0; j<all[dsItem(bb->genes,i)]->closure_rows; j++)
		{	
			if (!isregu[dsItem(all[dsItem(bb->genes,i)]->sequence,j)])
	                {
        	                isregu[dsItem(all[dsItem(bb->genes,i)]->sequence,j)] = TRUE;
                	        kk_1++;
                	}
			k = strlen(sequences[0])*dsItem(all[dsItem(bb->genes,i)]->sequence,j) + dsItem(all[dsItem(bb->genes,i)]->position,j);
			if (!isclo[k])
			{	
				if (k>0)
                                {
                                        isclo[k-1] = TRUE;
                                        if (k>1)  isclo[k-2] = TRUE;
                                }
                                if (k<s_rows*strlen(sequences[0])-1)
                                {
                                        isclo[k+1] = TRUE;
                                        if (k<s_rows*strlen(sequences[0])-2)  isclo[k+2] = TRUE;
                                }
				isclo[k] = TRUE;
				kk++;
			}
		}
	}
	for (i=0;i<s_rows*strlen(sequences[0]);i++)
		isclo[i] = FALSE;	
/**********************************************/
	fprintf (fw," Motif length: %d\n Motif number: %d\n Operon number: %d\n\n",all[dsItem(bb->genes,0)]->length,kk,kk_1);
	fprintf (fw,"------------------- Aligned Motif ------------------\n");
	if (po->ID) fprintf (fw,"#Seq\tposi\tID\tMotif\t\tAnnotation\n");
	else fprintf (fw,"#Motif\tSeq\tposi\tMotif\t\tAnnotation\n");
	sequences_regulon = alloc2c (kk,all[dsItem(bb->genes,0)]->length);
	int i3=0;
	int start=0/*, genome_num*/;
	int regulon_row=0,regulon_clo=0;
	for (i=0;i<(bb->block_rows);i++)
		i3 += get_num_TF (all[dsItem(bb->genes,i)]);
	char **reg_TF;
	reg_TF = alloc2c (i3,10);
	i3=0;
	for (i=0;i<(bb->block_rows);i++)
	{
		/*printf ("%d\t%d\n",dsItem(bb->genes,i),clo_matr1[dsItem(bb->genes,i)][0]);*/
/*		genome_num = get_genome_num_from_closure (all[dsItem(bb->genes,i)]);*/
		for (j=0; j<all[dsItem(bb->genes,i)]->closure_rows; j++)
		{
			k = strlen(sequences[0])*dsItem(all[dsItem(bb->genes,i)]->sequence,j) + dsItem(all[dsItem(bb->genes,i)]->position,j);
			if (!isclo[k])
			{
				/*make the scope (-2,2) adjacent to the current position TRUE*/
				isclo[k] = TRUE;
				if (k>0) 
				{
					isclo[k-1] = TRUE;
					if (k>1) isclo[k-2] = TRUE;
				}
				if (k<s_rows*strlen(sequences[0])-1)
				{
					isclo[k+1] = TRUE;
					if (k<s_rows*strlen(sequences[0])-2) isclo[k+2] = TRUE;
				}
				if (po->ID) fprintf (fw,">Motif-%d\t%d\t%d\t%d\t%s\t",num+1,dsItem(all[dsItem(bb->genes,i)]->sequence,j),dsItem(all[dsItem(bb->genes,i)]->position,j), dsItem(all[dsItem(bb->genes,i)]->position,j)+all[dsItem(bb->genes,0)]->length-1,locus_id[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]);
				else fprintf (fw,">Motif-%d\t%d\t%d\t%d\t",num+1,dsItem(all[dsItem(bb->genes,i)]->sequence,j),dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)], dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]+all[dsItem(bb->genes,0)]->length-1);
				/*adjsut closures base on p_opt and q_opt*/
/*				printf ("%d\t%d\t%d\t%d\t%d\n",dsItem(bb->genes,i),dsItem(bb->genes,0),clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)],dsItem(all[dsItem(bb->genes,i)]->sequence,j),dsItem(all[dsItem(bb->genes,i)]->position,j));*/
				if (dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)] < 0) 
				{
					start = 0;
					for (k=0;k<dsItem(all[dsItem(bb->genes,i)]->position,j)-clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)];k++)
					{
						fputc ('N',fw);
	                                        sequences_regulon[regulon_row][regulon_clo++] = 'N';
					}
					for (k=dsItem(all[dsItem(bb->genes,i)]->position,j)-clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)];k<all[dsItem(bb->genes,i)]->length;k++)
					{
						fprintf (fw,"%c",sequences[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k]);
	                                        sequences_regulon[regulon_row][regulon_clo++] = sequences[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k];
					}
				}
				else if ((dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]+all[dsItem(bb->genes,i)]->length) > strlen(sequences[0]))  
				{
					start =dsItem(all[dsItem(bb->genes,i)]->position,j);
					for (k=(dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]);k<strlen(sequences[0]);k++)
					{
						fprintf (fw,"%c",sequences[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k]);
	                                        sequences_regulon[regulon_row][regulon_clo++] = sequences[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k];
					}
					for (k=0;k<dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]+all[dsItem(bb->genes,i)]->length-strlen(sequences[0]);k++)
					{
						fputc ('N',fw);
                                                sequences_regulon[regulon_row][regulon_clo++]='N';
					}
				}
				else 
				{
					start = dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)];
					for (k=start; k< start+all[dsItem(bb->genes,i)]->length; k++)
					{
						fprintf (fw,"%c",sequences[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k]);
						sequences_regulon[regulon_row][regulon_clo++] = sequences[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k];
					}
				}
				regulon_row++;
				regulon_clo = 0;
				fprintf (fw, "\t");
				if (po->IS_SWITCH)
                	        {
                        	        int i1=0, i2=0;
                                	for (k = 0; k<anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->num;k++)
	                                {
        	                                i1 = dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->init,k);
                	                        i2 = dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->end,k);
                        	                if ((dsItem(all[dsItem(bb->genes,i)]->position,j)) >= i1&&(dsItem(all[dsItem(bb->genes,i)]->position,j)) <= i2)
                                	        {
                                        	        reg_TF[i3++] = anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->TF[k];
							fprintf (fw,"%s",anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->TF[k]);
                                                	fprintf (fw,"_%d",dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->init,k));
	                                                fprintf (fw,"_%d ",dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->end,k));
        	                                }
	
        	                        }
                	        }
				fprintf (fw, "\n");
			}
		}
	}
	if (po->IS_SWITCH) 
		print_frequency_anno (fw, reg_TF, i3,kk);
	if (po->IS_reference_H)
		print_operons (fw, sequences_regulon, genome,kk,oper_num_all);
	fprintf (fw,"----------------------------------------------------\n");
/*	for (i=0;i<ver;i++)
		free (arr_c1[i]);
	free (arr_c1);*/
}
/******************************************************************/

static void print_regulon_vertical( FILE *fw,  Block* bb, int num )
{
	int i,j,kk=0;
	long int k=0;
	verboseDot();
	fprintf (fw,"\n\n*********************************************************\n");
	fprintf (fw,"Candidate Regulon %d:\t ",(num+1));
	for (i=0;i< (bb->block_rows); i++)
		 fprintf (fw,"closure %d ",(dsItem(bb->genes,i)+1));
	fprintf (fw, "\n");
	fprintf (fw,"*********************************************************\n\n");
	
	for (i=0;i<(bb->block_rows);i++)
		kk += all[dsItem(bb->genes,i)]->closure_rows;
	fprintf (fw," Motif length: %d\n Motif number: %d\n",all[dsItem(bb->genes,0)]->length,kk);
	fprintf (fw,"------------------- Aligned Motif ------------------\n");
	if (po->ID) fprintf (fw,"#Genome\tSeq\tposi\tID\tMotif\t\tAnnotation\n");
	else fprintf (fw,"#Genome\tSeq\tposi\tMotif\t\tAnnotation\n");
	sequences_regulon = alloc2c (kk,all[dsItem(bb->genes,0)]->length);
	int i3=0;
	int start=0, genome_num;
	int regulon_row=0,regulon_clo=0;
	for (i=0;i<(bb->block_rows);i++)
		i3 += get_num_TF (all[dsItem(bb->genes,i)]);
	char **reg_TF;
	reg_TF = alloc2c (i3,10);
	i3=0;
	for (i=0;i<(bb->block_rows);i++)
	{
		/*printf ("%d\t%d\n",dsItem(bb->genes,i),clo_matr1[dsItem(bb->genes,i)][0]);*/
		genome_num = get_genome_num_from_closure (all[dsItem(bb->genes,i)]);
		for (j=0; j<all[dsItem(bb->genes,i)]->closure_rows; j++)
		{
		/*	if (po->ID) fprintf (fw,"%d\t%d\t%d\t%s\t",genome_num,dsItem(all[dsItem(bb->genes,i)]->sequence,j),dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)],locus_id[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]);
			else*/ fprintf (fw,"%d\t%d\t%d\t",genome_num,dsItem(all[dsItem(bb->genes,i)]->sequence,j),dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]);
			/*adjsut closures base on p_opt and q_opt*/
/*			printf ("%d\t%d\t%d\t%d\t%d\n",dsItem(bb->genes,i),dsItem(bb->genes,0),clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)],dsItem(all[dsItem(bb->genes,i)]->sequence,j),dsItem(all[dsItem(bb->genes,i)]->position,j));*/
			if (dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)] < 0) 
			{
				start = 0;
				for (k=0;k<dsItem(all[dsItem(bb->genes,i)]->position,j)-clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)];k++)
				{
					fputc ('N',fw);
                                        sequences_regulon[regulon_row][regulon_clo++] = 'N';
				}
				for (k=dsItem(all[dsItem(bb->genes,i)]->position,j)-clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)];k<all[dsItem(bb->genes,i)]->length;k++)
				{
					fprintf (fw,"%c",genome[genome_num]->sequences_r[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k]);
                                        sequences_regulon[regulon_row][regulon_clo++] = genome[genome_num]->sequences_r[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k];
				}
			}
			else if ((dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]+all[dsItem(bb->genes,i)]->length) > strlen(genome[genome_num]->sequences_r[0]))  
			{
				start =dsItem(all[dsItem(bb->genes,i)]->position,j);
				for (k=(dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]);k<strlen(genome[genome_num]->sequences_r[0]);k++)
				{
					fprintf (fw,"%c",genome[genome_num]->sequences_r[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k]);
                                        sequences_regulon[regulon_row][regulon_clo++] = genome[genome_num]->sequences_r[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k];
				}
				for (k=0;k<dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)]+all[dsItem(bb->genes,i)]->length-strlen(genome[genome_num]->sequences_r[0]);k++)
				{
					fputc ('N',fw);
                                        sequences_regulon[regulon_row][regulon_clo++]='N';
				}
			}
			else 
			{
				start = dsItem(all[dsItem(bb->genes,i)]->position,j)+clo_matr1[dsItem(bb->genes,i)][dsItem(bb->genes,0)];
				for (k=start; k< start+all[dsItem(bb->genes,i)]->length; k++)
				{
					fprintf (fw,"%c",genome[genome_num]->sequences_r[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k]);
					sequences_regulon[regulon_row][regulon_clo++] = genome[genome_num]->sequences_r[dsItem(all[dsItem(bb->genes,i)]->sequence,j)][k];
				}
			}
			regulon_row++;
			regulon_clo = 0;
			fprintf (fw, "\t");
			if (po->IS_SWITCH)
               	        {
                       	        int i1=0, i2=0;
                               	for (k = 0; k<anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->num;k++)
                                {
       	                                i1 = dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->init,k);
               	                        i2 = dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->end,k);
                       	                if ((dsItem(all[dsItem(bb->genes,i)]->position,j)) >= i1&&(dsItem(all[dsItem(bb->genes,i)]->position,j)) <= i2)
                               	        {
                                       	        reg_TF[i3++] = anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->TF[k];
						fprintf (fw,"%s",anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->TF[k]);
                                               	fprintf (fw,"_%d",dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->init,k));
                                                fprintf (fw,"_%d ",dsItem(anno[dsItem(all[dsItem(bb->genes,i)]->sequence,j)]->end,k));
       	                                }	
       	                        }
               	        }
			fprintf (fw, "\n");
		}
	}
	if (po->IS_SWITCH) 
		print_frequency_anno (fw, reg_TF, i3,kk);
	if (po->IS_reference)
		print_operons (fw, sequences_regulon, genome,kk,oper_num_all);	
	fprintf (fw,"----------------------------------------------------\n");
/*	for (i=0;i<ver;i++)
		free (arr_c1[i]);
	free (arr_c1);*/
}
/******************************************************************/
