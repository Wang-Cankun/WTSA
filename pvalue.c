/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 16, 2010
 * all the functions realted to calculate pvalue
 */
#include "pvalue.h"

/***********************************************************************************/
continuous **markov (char **sequences_r, int seq_num)
{
	size_t i,j;
	continuous **markov;
	markov= alloc2dd (5,5);
        for (i=0;i<5;i++)
                for (j=0;j<5;j++)
                        markov[i][j]=0;
        size_t flag_markov=0;
        sum_markov=0;
	size_t length_temp = strlen(sequences_r[0]);
        for (i=0;i<seq_num;i++)
	{
		/*strlen() return a size_t format so if we compare it with int (for (j=1;j<strlen(sequences_r[i]);j++)), there will be segemental fault
		 * here we set strlen(sequences_r[0]) to a int and do the following steps*/
		length_temp = strlen(sequences_r[i]);
		for (j=1;j<length_temp;j++)
                {
                        sum_markov++;
                        if (sequences_r[i][j-1]=='A' || sequences_r[i][j-1]=='a') flag_markov=1;
                        else if (sequences_r[i][j-1]=='G' || sequences_r[i][j-1]=='g') flag_markov=2;
                        else if (sequences_r[i][j-1]=='C' || sequences_r[i][j-1]=='c') flag_markov=3;
                        else if (sequences_r[i][j-1]=='T' || sequences_r[i][j-1]=='t') flag_markov=4;
                        markov[flag_markov][0]++;
			if (sequences_r[i][j]=='A' || sequences_r[i][j]=='a') markov[flag_markov][1]++;
                        else if (sequences_r[i][j]=='G'||sequences_r[i][j]=='g') markov[flag_markov][2]++;
                        else if (sequences_r[i][j]=='C'||sequences_r[i][j]=='c') markov[flag_markov][3]++;
                        else if (sequences_r[i][j]=='T'||sequences_r[i][j]=='t') markov[flag_markov][4]++;
                }
	}
        for (i=1;i<5;i++)
                for (j=1;j<5;j++)
                        markov[i][j]=markov[i][j]/markov[i][0];
        for (i=0;i<seq_num;i++)
        {
                length_temp = strlen(sequences_r[i]);
		if (sequences_r[i][length_temp-1]=='A' || sequences_r[i][length_temp-1]=='a') flag_markov=1;
                else if (sequences_r[i][length_temp-1]=='G' || sequences_r[i][length_temp-1]=='g') flag_markov=2;
                else if (sequences_r[i][length_temp-1]=='C' || sequences_r[i][length_temp-1]=='c') flag_markov=3;
                else if (sequences_r[i][length_temp-1]=='T' || sequences_r[i][length_temp-1]=='t') flag_markov=4;
		markov[flag_markov][0]++;
                sum_markov++;
        }
        for (i=1;i<5;i++)
                markov[i][0]=markov[i][0]/sum_markov;


	/*for (i=0;i<5;i++){
                for (j=0;j<5;j++){
                        printf("%f\t",markov[i][j]);
		}
		printf("\n");
	}*/
	return markov;
}



/*WTSA 3-order markov matrix*/
continuous **d_markov (char **sequences_r, int seq_num)
{
	size_t i,j;
	continuous **d_markov;
	d_markov= alloc2dd (17,5);
        for (i=0;i<17;i++)
                for (j=0;j<5;j++)
                        d_markov[i][j]=0;
        size_t flag_markov=0;
        sum_markov=0;
	size_t length_temp = strlen(sequences_r[0]);
        for (i=0;i<seq_num;i++)
	{
		/*strlen() return a size_t format so if we compare it with int (for (j=1;j<strlen(sequences_r[i]);j++)), there will be segemental fault
		 * here we set strlen(sequences_r[0]) to a int and do the following steps*/
		length_temp = strlen(sequences_r[i]);
		for (j=2;j<length_temp;j++)
                {
                        sum_markov++;
                        if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]=='A') flag_markov=1;
                        else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='a') flag_markov=1;
			else if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]=='T') flag_markov=2;
                        else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='t') flag_markov=2;
			else if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]=='C') flag_markov=3;
                        else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='c') flag_markov=3;
			else if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]=='G') flag_markov=4;
                        else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='g') flag_markov=4;

                        else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='A') flag_markov=5;
                        else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='a') flag_markov=5;
			else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='T') flag_markov=6;
                        else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='t') flag_markov=6;
			else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='C') flag_markov=7;
                        else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='c') flag_markov=7;
			else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='G') flag_markov=8;
                        else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='g') flag_markov=8;

                        else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='A') flag_markov=9;
                        else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='a') flag_markov=9;
			else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='T') flag_markov=10;
                        else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='t') flag_markov=10;
			else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='C') flag_markov=11;
                        else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='c') flag_markov=11;
			else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='G') flag_markov=12;
                        else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='g') flag_markov=12;

                        else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='A') flag_markov=13;
                        else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='a') flag_markov=13;
			else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='T') flag_markov=14;
                        else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='t') flag_markov=14;
			else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='C') flag_markov=15;
                        else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='c') flag_markov=15;
			else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='G') flag_markov=16;
                        else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='g') flag_markov=16;


                        d_markov[flag_markov][0]++;
			if (sequences_r[i][j]=='A' || sequences_r[i][j]=='a') d_markov[flag_markov][1]++;
                        else if (sequences_r[i][j]=='G'||sequences_r[i][j]=='g') d_markov[flag_markov][2]++;
                        else if (sequences_r[i][j]=='C'||sequences_r[i][j]=='c') d_markov[flag_markov][3]++;
                        else if (sequences_r[i][j]=='T'||sequences_r[i][j]=='t') d_markov[flag_markov][4]++;
                }
	}
	
        for (i=1;i<17;i++)
                for (j=1;j<5;j++)
                        d_markov[i][j]=d_markov[i][j]/d_markov[i][0];
        for (i=0;i<seq_num;i++)
        {
                j = strlen(sequences_r[i]);
		

	if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]== 'A') flag_markov=1;
else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='a') flag_markov=1;
else if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]=='G') flag_markov=2;
else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='g') flag_markov=2;
else if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]=='C') flag_markov=3;
else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='c') flag_markov=3;
else if (sequences_r[i][j-2]=='A' && sequences_r[i][j-1]=='T') flag_markov=4;
else if (sequences_r[i][j-2]=='a' && sequences_r[i][j-1]=='t') flag_markov=4;

else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='A') flag_markov=5;
else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='a') flag_markov=5;
else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='G') flag_markov=6;
else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='g') flag_markov=6;
else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='C') flag_markov=7;
else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='c') flag_markov=7;
else if (sequences_r[i][j-2]=='G' && sequences_r[i][j-1]=='T') flag_markov=8;
else if (sequences_r[i][j-2]=='g' && sequences_r[i][j-1]=='t') flag_markov=8;

else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='A') flag_markov=9;
else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='a') flag_markov=9;
else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='G') flag_markov=10;
else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='g') flag_markov=10;
else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='C') flag_markov=11;
else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='c') flag_markov=11;
else if (sequences_r[i][j-2]=='C' && sequences_r[i][j-1]=='T') flag_markov=12;
else if (sequences_r[i][j-2]=='c' && sequences_r[i][j-1]=='t') flag_markov=12;


else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='A') flag_markov=13;
else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='a') flag_markov=13;
else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='G') flag_markov=14;
else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='g') flag_markov=14;
else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='C') flag_markov=15;
else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='c') flag_markov=15;
else if (sequences_r[i][j-2]=='T' && sequences_r[i][j-1]=='T') flag_markov=16;
else if (sequences_r[i][j-2]=='t' && sequences_r[i][j-1]=='t') flag_markov=16;


		d_markov[flag_markov][0]++;
                sum_markov++;
        }
        for (i=1;i<17;i++)
                d_markov[i][0]=d_markov[i][0]/sum_markov;


	/*for (i=0;i<17;i++){
                for (j=0;j<5;j++){
                        printf("%f\t",d_markov[i][j]);
		}
		printf("\n");
	}*/
	return d_markov;
}

/***********************************************************************************/
continuous aver_score_closure(char **sequences_2, continuous **scoreM, continuous *score, int motif_number, int length_local_1)
{
	int i,j;
	continuous AveScore=0;
	for (i=0;i<motif_number;i++)
        {
        	score[i]=0;
                for (j=0;j<length_local_1;j++)
                {
                	if (sequences_2[i][j]=='A' || sequences_2[i][j]=='a') score[i]=score[i]+scoreM[1][j];
                        else if (sequences_2[i][j]=='G' || sequences_2[i][j]=='g') score[i]=score[i]+scoreM[2][j];
                        else if (sequences_2[i][j]=='C' || sequences_2[i][j]=='c') score[i]=score[i]+scoreM[3][j];
                        else if (sequences_2[i][j]=='T' || sequences_2[i][j]=='t') score[i]=score[i]+scoreM[4][j];
                }
	        AveScore = AveScore + score[i];
        }
        AveScore = AveScore / motif_number;		
	return AveScore;
}	
/***********************************************************************************/
int motif_num_closure (Reference_genome *genome_1, int motif_length, continuous **scoreM, continuous thre)
{
	int i,j,k;
	continuous score_scan=0;
	int motif_number=0;
        for (i=0;i<genome_1->seq_num;i++)
	{	
		for (j=0;j< (genome_1->clo_num[i])-motif_length+1;j++)
                {
	                score_scan=0;
                        for (k=0;k<motif_length;k++)
                        {
                       		if (genome_1->sequences_r[i][k+j]=='A' || genome_1->sequences_r[i][k+j]=='a')
                                	score_scan = score_scan + scoreM[1][k];
                                if (genome_1->sequences_r[i][k+j]=='G' || genome_1->sequences_r[i][k+j]=='g')
                                	score_scan = score_scan + scoreM[2][k];
                                if (genome_1->sequences_r[i][k+j]=='C' || genome_1->sequences_r[i][k+j]=='c')
                                	score_scan = score_scan + scoreM[3][k];
                                if (genome_1->sequences_r[i][k+j]=='T' || genome_1->sequences_r[i][k+j]=='t')
                                	score_scan = score_scan + scoreM[4][k];
                        }
                        if (score_scan>thre && motif_number < CLOSURE_SIZE-10)
                                        motif_number++;
		}
	}
	return motif_number;
}
/***********************************************************************************/
continuous *motif_num_R_closure (Reference_genome *genome_1, int motif_length, continuous **scoreM, continuous thre)
{
	int i,j,k,length_ave=0,numberfore=1,randomNum,simulation=po->simu,Rt,num_all_R=0;
        continuous length_ave_1=0, score_scan;
	continuous pp[5];

	continuous  motif_num_R_temp=0;
	continuous *motif_num_R;
	AllocArray (motif_num_R,2);
	motif_num_R[0] = motif_num_R[1] =0;

	for (i=0;i<genome_1->seq_num;i++)
                length_ave_1=length_ave_1+genome_1->clo_num[i];
        length_ave=ceil(length_ave_1/genome_1->seq_num);
        srand((unsigned)time(NULL));
        char *randomdata;
        AllocArray (randomdata,length_ave);
        for (Rt=0;Rt<30*simulation;Rt++)
	{
                for (i=0;i<4;i++)
                {
                        for(j=0;j<length_ave;j++)
                        {
                                num_all_R++;
                                for (k=1;k<5;k++)
                                        pp[k]=genome_1->markov_matrix[numberfore][k];
                                randomNum=rand()%100;
                                if (randomNum<pp[1]*100) {randomdata[j]='A';numberfore=1;}
                                else if (randomNum<(pp[1]+pp[2])*100) {randomdata[j]='G';numberfore=2;}
                                else if (randomNum<(pp[1]+pp[2]+pp[3])*100) {randomdata[j]='C';numberfore=3;}
                                else {randomdata[j]='T';numberfore=4;}
                        }
                        for(j=0;j<length_ave-motif_length+1;j++)
                        {
                                score_scan=0;
                                for (k=0;k<motif_length;k++)
                                {
                                        if (randomdata[k+j]=='A' || randomdata[k+j]=='a')
                                                score_scan=score_scan+scoreM[1][k];
                                        else if (randomdata[k+j]=='G' || randomdata[k+j]=='g')
                                                score_scan=score_scan+scoreM[2][k];
                                        else if (randomdata[k+j]=='C' || randomdata[k+j]=='c')
                                                score_scan=score_scan+scoreM[3][k];
                                        else if (randomdata[k+j]=='T' || randomdata[k+j]=='t')
                                                score_scan=score_scan+scoreM[4][k];
                                }
                                if (score_scan>thre)
                                	motif_num_R_temp++;
                        }
                }
	}
	if (motif_num_R_temp == 0) motif_num_R_temp++;
/*	printf ("*****%3.2f\n",(4*(Rt))/(motif_num_R_temp*genome_1->seq_num));*/
	motif_num_R[0] = (motif_num_R_temp*genome_1->seq_num)/(4*(Rt));
	motif_num_R[1] = (4*(Rt))/(motif_num_R_temp*genome_1->seq_num);
/*	printf ("%d\t%d\n",motif_num_R[0],motif_num_R[1]);*/
	/*free(randomdata);*/
	return motif_num_R;
}

/***********************************************************************************/
continuous get_pvalue(Reference_genome *genome_1, int motif_num, continuous *motif_num_R)
{
/*	printf ("%d\t%3.2f\n",motif_num,motif_num_R[0]);*/
	int i,j;
	continuous poisson,one_tmp=1,pvalue;
	pvalue=0;
        while (motif_num_R[0] >600)
	{
		motif_num_R[0] = motif_num_R[0]/2;
		motif_num = motif_num/2;
	}
	poisson=one_tmp/exp(motif_num_R[0]);
	/*if po->approximate j=1*/
	if (po->approximate) j=1;
	else j=300;
        for (i=0;i<(int)(motif_num)+j;i++)
        {
        	if (i>(int)(motif_num)-1) pvalue=pvalue+poisson;
                poisson=poisson*motif_num_R[0]/(i+1);
	}
	return pvalue;
}
/***********************************************************************************/
void print_operons (FILE *fw, char **sequences_regulon, Reference_genome **genome, int motif_num, int oper_num)
{
	int i,j;
	int motif_number=0;
       	continuous *motif_number_R;
	continuous **scoreR, pvalue=1, pvalue_temp, ave_socre,*score, thre;
	AllocArray (score, motif_num);
	double length_local = po->MOTIFLENGTH;
        
	/*get_profile function*/
	discrete **frequency;
        scoreR = alloc2dd (5,po->MOTIFLENGTH);
        for (i=0; i<5; i++)
                for (j=0; j<po->MOTIFLENGTH; j++)
                        scoreR[i][j]=0;
        frequency = frequency_matrix (sequences_regulon, 5, po->MOTIFLENGTH, motif_num);
        for (i=1;i<5;i++)
                for (j=0;j<po->MOTIFLENGTH;j++)
                        scoreR[i][j]= log((frequency[i][j]+p_markov[i][0])/(p_markov[i][0]*(motif_num+1)));
	/*improve the profile*/	
	if (po->palindromic)    scoreR = impovre_profle_palindromic (scoreR, length_local, sequences_regulon, 5 ,motif_num);
	else scoreR = impovre_profle (scoreR, length_local);

	/*get the average score*/
	ave_socre = aver_score_closure(sequences_regulon, scoreR, score, motif_num, po->MOTIFLENGTH);

	fprintf (fw,"----------------------------------------------------\nPredicted_operons:\t");
	for (i=0;i<oper_num_all; i++)
	{
	/*	if ((i+1)%10==0) verboseDot();*/
		genome[i]->markov_matrix = alloc2dd (5,5);
		genome[i]->markov_matrix = markov (genome[i]->sequences_r,genome[i]->seq_num);
		/* select the minimal pvalue for 0.3-0.9*/
		for (thre=0.3; thre<1;thre=thre+0.1)
		{
			motif_number = motif_num_closure (genome[i], po->MOTIFLENGTH , scoreR, thre*ave_socre);
			motif_number_R = motif_num_R_closure(genome[i], po->MOTIFLENGTH, scoreR, thre*ave_socre);
/*			printf ("%d\t%d\n",motif_number , motif_number_R);*/
			pvalue_temp = get_pvalue(genome[i], motif_number, motif_number_R);
/*			printf ("%d\t%s\t%d\t%d\t%3.2e\n",i,genome[i]->oper_name,motif_number,motif_number_R[0],pvalue_temp);*/
			if (pvalue_temp < pvalue) pvalue = pvalue_temp;
		}
	/*	printf ("***********%s\t%3.2e\n",genome[i]->oper_name,pvalue);*/
		genome[i]->score = -(100*log(pvalue));
		genome[i]->pvalue = pvalue;
/*		printf ("%3.2e\t",genome[i]->pvalue);
		if (pvalue < po->thre_pvalue) fprintf (fw, "%s(%3.2e)\t",genome[i]->oper_name,pvalue);*/
		pvalue_temp = pvalue = 1;
	}
	sort_genomes_list(genome, oper_num_all);	
/*	printf ("%3.2e\t%3.2e\n",genome[0]->pvalue,genome[1]->pvalue);*/
	int min = MIN ( motif_num, oper_num_all);
	for (i=0;i< min; i++)
	{
		if (genome[i]->pvalue < po->thre_pvalue)
			fprintf (fw, "%s(%3.2e)\t",genome[i]->oper_name,genome[i]->pvalue);
	}
	fprintf (fw,"\n");
/*	printf ("\n");*/
	free (score);
}
/***********************************************************************************/
int block_cmpr_g(const void *a, const void *b)
/*compare function for qsort, decreasing by score*/
{
        return ((*(Reference_genome **)b)->score - (*(Reference_genome **)a)->score);
}

void sort_genomes_list(Reference_genome **el, int n)
{
        qsort(el, n, sizeof *el, block_cmpr_g);
}

/************************************************************************************/
bool is_two_closures_regulon (Closures *a, Closures *b, Reference_genome *aa, Reference_genome *bb, int num_a, int num_b)
{
/*	printf ("1111%d\t%d\n",dsItem(all[0]->sequence,0),dsItem(all[0]->position,0));*/
	int i,k,motif_number=0;
	size_t j;
	char **sequence_combine;
	continuous **profile, thre, *motif_number_R, *motif_number_R_aa, *motif_number_R_bb, *score, pvalue_temp, pvalue =1,ave_socre;
	sequence_combine = alloc2c((a->closure_rows+b->closure_rows),po->MOTIFLENGTH);
	AllocArray (motif_number_R,2);
	AllocArray (score, a->closure_rows+b->closure_rows);
/*	printf ("%d\t%d\t%d\t%d\n",num_a,num_b,aa->seq_num,bb->seq_num);*/
	/*get the sequence file from the two closures*/
	for (i=0; i<a->closure_rows;i++)
	{
	/*	printf ("%d\t%d\n",dsItem(a->sequence,i),dsItem(a->position,i));*/
		k = dsItem(a->position,i);
		for (j=0;j<a->length+extend_len;j++)
		{
			sequence_combine[i][j] = aa->sequences_r[dsItem(a->sequence,i)][k++];
		}	
	}
/*	printf ("2222%d\t%d\n",dsItem(all[0]->sequence,0),dsItem(all[0]->position,0));*/
	for (i=0;i<b->closure_rows;i++)
	{
		/*printf ("%d\t%d\t%d\t%c\n",b->closure_rows,dsItem(b->position,i),clo_matr1[num_b][num_a],bb->sequences_r[0][0]);*/
		if (dsItem(b->position,i)+clo_matr1[num_b][num_a] < 0)
                {
			k = dsItem(b->position,i);
                        for (j=0;j<dsItem(b->position,i)-clo_matr1[num_b][num_a];j++)
				sequence_combine[i+a->closure_rows][j] = 'N';
                        for (j=dsItem(b->position,i)-clo_matr1[num_b][num_a];j<b->length+extend_len;j++)
                                sequence_combine[i+a->closure_rows][j] = bb->sequences_r[dsItem(b->sequence,i)][k++];                                
		}
		else if (dsItem(b->position,i)+clo_matr1[num_b][num_a]+b->length+extend_len >strlen(bb->sequences_r[0]))
                {
			k = dsItem(b->position,i);
                        for (j=0; j<strlen(bb->sequences_r[0])-dsItem(b->position,i)-clo_matr1[num_b][num_a];j++)
                                sequence_combine[i+a->closure_rows][j] = bb->sequences_r[dsItem(b->sequence,i)][k++];                            
			for (j=strlen(bb->sequences_r[0])-dsItem(b->position,i)-clo_matr1[num_b][num_a];j<b->length+extend_len;j++)
				sequence_combine[i+a->closure_rows][j] = 'N';
		}
                else 
                {
/*	                printf ("11111111111\t%d\t%d\n",dsItem(b->sequence,i),dsItem(b->position,i));*/
			k = dsItem(b->position,i);
			for (j=0;j<b->length+extend_len;j++)
			{
                                /*printf ("%c\n",bb->sequences_r[dsItem(b->sequence,i)][dsItem(b->position,i)]);*/
				sequence_combine[i+a->closure_rows][j] = bb->sequences_r[dsItem(b->sequence,i)][k++];   
			}
		}
	}
/*	printf ("3333%d\t%d\n",dsItem(all[0]->sequence,0),dsItem(all[0]->position,0));*/
	/*get the profile of sequence_combine*/	
	profile = get_profile (sequence_combine, 5, a->length+extend_len, a->closure_rows+b->closure_rows);
	/*improve the profile*/	
	if (po->palindromic)    profile = impovre_profle_palindromic (profile, b->length+extend_len, sequence_combine, 5 ,a->closure_rows+b->closure_rows);
	else  profile = impovre_profle (profile, b->length+extend_len);
	/*get the average score*/
        ave_socre = aver_score_closure(sequence_combine, profile, score, a->closure_rows+b->closure_rows, po->MOTIFLENGTH);
	/*calculate the pvalue of combined clousres*/
	aa->markov_matrix = alloc2dd (5,5);
        aa->markov_matrix = markov (aa->sequences_r,aa->seq_num);
	bb->markov_matrix = alloc2dd (5,5);
        bb->markov_matrix = markov (bb->sequences_r,bb->seq_num);
        for (thre=0.3; thre<1;thre=thre+0.1)
        {
		motif_number = motif_number_R[0]= motif_number_R[1]=0;
		motif_number += motif_num_closure (aa, po->MOTIFLENGTH , profile, thre*ave_socre);
		motif_number += motif_num_closure (bb, po->MOTIFLENGTH , profile, thre*ave_socre);
                motif_number_R_aa = motif_num_R_closure(aa, po->MOTIFLENGTH, profile, thre*ave_socre);
		motif_number_R_bb = motif_num_R_closure(bb, po->MOTIFLENGTH, profile, thre*ave_socre);
		motif_number_R[0] = motif_number_R_aa[0]+ motif_number_R_bb[0];
		motif_number_R[1] = motif_number_R_aa[1]+ motif_number_R_bb[1];
                pvalue_temp = get_pvalue(aa, motif_number, motif_number_R);
                if (pvalue_temp < pvalue) pvalue = pvalue_temp;
        }
/*	printf ("%3.2e\n",pvalue);*/
	
	/*free the memory*/
	free (motif_number_R);
	free (score);
	for (i=0;i<5;i++)
		free (profile[i]);
	free (profile);
/*	printf ("%d\t%d\t%3.2e\t%3.2e\t%3.2e\n",num_a,num_b,pvalue, a->significance, b->significance);*/
	if (pvalue <= MAX(a->significance,b->significance)) return TRUE;
	else return FALSE;
}

/************************************************************************************/
