/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 16, 2010     
 * the main process of the input fasta data: three times compare sequences
 */

#include "compare_sequence.h"
#define INF 5000000.00

/***********************************************************************/
/* set the initial value for adjent matrix generated by compare_sequences */
void init_dis()
{
	ver = po->TOPVERTICES*s_rows;
	AllocArray(genes, ver);
	AllocArray(conds, ver);
	int i;
	for (i=0;i<ver;i++)
	{
		genes[i]=i; 
		conds[i]=i;
	}
	rows=cols=ver; 
}
/***********************************************************************/
/* set the initial value for adjent matrix generated by get_closure_matrix_1 if po->IS_closure is TRUE*/
void init_dis_1()
{
	ver = clo_num;
	AllocArray(genes, ver);
	AllocArray(conds, ver);
	int i;
	for (i=0;i<ver;i++)
	{
		genes[i]=i; 
		conds[i]=i;
	}
	rows=cols=ver; 
}
/***********************************************************************/
/* some default setting for input length 5-20 */
void parameter_choice()
{
	if (po->MOTIFLENGTH == 5  && po->end_weight == 2) {po->threshold_1=3;po->threshold_2=4;po->threshold=2;}
	if (po->MOTIFLENGTH ==5  && po->end_weight == 1.5) {po->threshold_1=2;po->threshold_2=3;po->threshold=2;}
	if (po->MOTIFLENGTH ==6  && po->end_weight == 2) {po->threshold_1=3;po->threshold_2=4;po->threshold=2;}
	if (po->MOTIFLENGTH ==6  && po->end_weight == 1.5) {po->threshold_1=2;po->threshold_2=3;po->threshold=2;}
	if (po->MOTIFLENGTH ==7  && po->end_weight == 2) {po->threshold_1=4;po->threshold_2=5;po->threshold=3;}
	if (po->MOTIFLENGTH ==7  && po->end_weight == 1.5) {po->threshold_1=3;po->threshold_2=4;po->threshold=3;}
	if (po->MOTIFLENGTH ==8  && po->end_weight == 2) {po->threshold_1=4;po->threshold_2=5;po->threshold=3;}
	if (po->MOTIFLENGTH ==8  && po->end_weight == 1.5) {po->threshold_1=3;po->threshold_2=4;po->threshold=3;}
	if (po->MOTIFLENGTH ==9  && po->end_weight == 2) {po->threshold_1=5;po->threshold_2=6;po->threshold=4;}
	if (po->MOTIFLENGTH ==9  && po->end_weight == 1.5) {po->threshold_1=4;po->threshold_2=5;po->threshold=4;}
	if (po->MOTIFLENGTH ==10  && po->end_weight == 2) {po->threshold_1=11;po->threshold_2=14;po->threshold=5;}
	if (po->MOTIFLENGTH ==10  && po->end_weight == 1.5) {po->threshold_1=6;po->threshold_2=8;po->threshold=5;}
	if (po->MOTIFLENGTH ==11  && po->end_weight == 2) {po->threshold_1=11;po->threshold_2=14;po->threshold=6;}
	if (po->MOTIFLENGTH ==11  && po->end_weight == 1.5) {po->threshold_1=7;po->threshold_2=10;po->threshold=6;}
	if (po->MOTIFLENGTH ==12  && po->end_weight == 2) {po->threshold_1=12;po->threshold_2=14;po->threshold=6;}
	if (po->MOTIFLENGTH ==12  && po->end_weight == 1.5) {po->threshold_1=9;po->threshold_2=12;po->threshold=6;}
	if (po->MOTIFLENGTH ==13  && po->end_weight == 2) {po->threshold_1=13;po->threshold_2=16;po->threshold=7;}
	if (po->MOTIFLENGTH ==13  && po->end_weight == 1.5) {po->threshold_1=10;po->threshold_2=13;po->threshold=7;}
	if (po->MOTIFLENGTH ==14  && po->end_weight == 2) {po->threshold_1=15;po->threshold_2=17;po->threshold=7;}
	if (po->MOTIFLENGTH ==14  && po->end_weight == 1.5) {po->threshold_1=10;po->threshold_2=14;po->threshold=7;}
	if (po->MOTIFLENGTH ==14  && po->end_weight == 1) {po->threshold_1=10;po->threshold_2=10;po->threshold=7;}
	if (po->MOTIFLENGTH ==15  && po->end_weight == 2) {po->threshold_1=16;po->threshold_2=19;po->threshold=7;}
	if (po->MOTIFLENGTH ==15  && po->end_weight == 1.5) {po->threshold_1=14;po->threshold_2=15;po->threshold=7;}
	if (po->MOTIFLENGTH ==16  && po->end_weight == 2) {po->threshold_1=16;po->threshold_2=19;po->threshold=8;}
	if (po->MOTIFLENGTH ==16  && po->end_weight == 1.5) {po->threshold_1=14;po->threshold_2=16;po->threshold=8;}
	if (po->MOTIFLENGTH ==17  && po->end_weight == 2) {po->threshold_1=16;po->threshold_2=18;po->threshold=8;}
	if (po->MOTIFLENGTH ==17  && po->end_weight == 1.5) {po->threshold_1=13;po->threshold_2=15;po->threshold=8;}
	if (po->MOTIFLENGTH ==18  && po->end_weight == 2) {po->threshold_1=16;po->threshold_2=19;po->threshold=9;}
	if (po->MOTIFLENGTH ==18  && po->end_weight == 1.5) {po->threshold_1=14;po->threshold_2=16;po->threshold=9;}
	if (po->MOTIFLENGTH ==19  && po->end_weight == 2) {po->threshold_1=16;po->threshold_2=18;po->threshold=9;}
	if (po->MOTIFLENGTH ==19  && po->end_weight == 1.5) {po->threshold_1=13;po->threshold_2=15;po->threshold=9;}
	if (po->MOTIFLENGTH ==20  && po->end_weight == 2) {po->threshold_1=16;po->threshold_2=18;po->threshold=9;}
	if (po->MOTIFLENGTH ==20  && po->end_weight == 1.5) {po->threshold_1=13;po->threshold_2=15;po->threshold=9;}
}
/***********************************************************************/
static bool **matrix_no_continuous_equal(bool **matrix)
{
	/* consider 4 continuous equal subsequence*/
	int i,j,pos;
	for(i=0;i<s_rows;i++)
	{
	        for (j=0;j<s_cols;j++)
                {
                        if (sequences[i][j]==sequences[i][j+1] && sequences[i][j+1]==sequences[i][j+2] && sequences[i][j+2]==sequences[i][j+3])
                        {
                                matrix[i][j]=TRUE;
                                continue;
                        }
                        if (sequences[i][j+po->MOTIFLENGTH-4]==sequences[i][j+po->MOTIFLENGTH-3] && sequences[i][j+po->MOTIFLENGTH-3]==sequences[i][j+po->MOTIFLENGTH-2] && sequences[i][j+po->MOTIFLENGTH-2]==sequences[i][j+po->MOTIFLENGTH-1])
                        {
                                matrix[i][j]=TRUE;
                                continue;
                        }
                        for (pos=0;pos<po->MOTIFLENGTH-4;pos++)
			{
                                if (sequences[i][j+pos]==sequences[i][j+pos+1] && sequences[i][j+pos+1]==sequences[i][j+pos+2] && sequences[i][j+pos+2]==sequences[i][j+pos+3] && sequences[i][j+pos+3]==sequences[i][j+pos+4])
                                {
                                        matrix[i][j]=TRUE;
                                        break;
                                }
			}
                }
	}
	return matrix;
}
/***********************************************************************/
static void pairwise_comparison_first ( bool **matrix3, int **matrix, bool *match, bool *match1, int lower, int upper, struct dispos *marray, struct edgepos *marray1, int mpos, int msize, int maxscore, int maxindex, int *startpos)
{
	
	int i,j,p,q,k;
	int temp_largest;
	int max1;
	int max2;
	float max_all[s_rows];
	uglyTime("compare first start", i);
	continuous **d1, **d2;
        d1 = alloc2dd (s_cols,s_cols);
        d2 = alloc2dd (s_cols,s_cols);

	maxindex = 0; 
	if (!po->middle_enhance && !po->no_enhance)
	printf ("\nMotif length is: %d\nSignal strengthening started\n", po->MOTIFLENGTH);
	fflush(stdout);
	/*printf("1%d\t\n",s_rows);*/
	
	for(i=0;i<s_rows;i++)
	{
		max_all[i] = height_matrix[i][0];
		for (temp_largest = 1; temp_largest < s_col[i]; temp_largest++)
		{
			
			if (height_matrix[i][temp_largest] > max_all[i] )
			{
			max_all[i] = height_matrix[i][temp_largest];
			max_all[i] = logf((float)max_all[i]);
			}
		}
		/*printf("%d\t%f\t\n",i,max_all[i] );*/
	}

	for(i=0;i<s_rows;i++)
	{
		/* one dot represent ten sequences */
		if ((i+1)%10==0) verboseDot(); 
		/*uglyTime("first for i=%d", i);*/

		for(p=i;p<s_rows;p++)
		{		
	/*		max2 = height_matrix[p][0];
			max1 = height_matrix[i][0];
		for (temp_largest = 1; temp_largest < s_col[i]; temp_largest++)
		{
			printf("%d\t%d\t%d\t%d\t%d\t\n",i,j,p,q,max1);
			if (height_matrix[i][temp_largest] > max1)
			{
			max1 = height_matrix[i][temp_largest];
			}
				
		}
			for (temp_largest = 1; temp_largest < s_col[p]; temp_largest++)
			{
				if (height_matrix[p][temp_largest] > max2)
				{
				max2 = height_matrix[p][temp_largest];
				}
					
			}*/
			/*initialize the values of d1*/
			for (j=0;j<s_col[i];j++)
				for (q=0;q<s_col[p];q++)
					{ d1[j][q]=0; d2[j][q]=0;}
			/*store values of d1 for the ith sequence of first alignment*/
			for(j=0;j<s_col[i];j++)
			{
				if (matrix3[i][j])   continue; /*check matrix_no_continuous_equal*/
				for(q=0;q<s_col[p];q++)
				{

					if (matrix3[p][q]) continue; /*check matrix_no_continuous_equal*/
					/*do not compare with itself*/
					if (q==j && i==p) { d1[j][q]=0; continue;}
					
					/*if (j>0 && q>0 ) d1[j][q] = d1[j-1][q-1] - fre_matrix[p][seq_matrix[i][j-1]][q-1] + fre_matrix[p][seq_matrix[i][j+po->MOTIFLENGTH-1]][q+po->MOTIFLENGTH-1];
					else*/ d1[j][q]=get_similarity_between_two_patterns (i,p,j,q,po->MOTIFLENGTH);
					
					/* set a threshold for d1 */
					if (d1[j][q] ==0.0) {continue;}
					if (max2==1) {d2[j][q]=0.0;}
					else {
						d2[j][q] = improve_similarity_between_two_patterns(i, p, j, q, lower, upper, d1[j][q])* (log1pf((float)height_matrix[i][j])/max_all[i]+log1pf((float)height_matrix[p][q])/max_all[p]);
						/*d2[j][q] = improve_similarity_between_two_patterns(i, p, j, q, lower, upper, d1[j][q])* (log1pf(height_matrix[i][j])/logf(max1)+log1pf(height_matrix[p][q])/logf(max2));*/
						/*d2[j][q] = improve_similarity_between_two_patterns(i, p, j, q, lower, upper, d1[j][q]);*/
						
						}
        			}
					
      			}

			/* initialize match and match1 */
			for (j=0;j<s_cols;j++) { match[j]=FALSE; match1[j]=FALSE;}
			msize = 2*(s_col[i]+s_col[p]-1)-2;
			/*select top2 in each alignment of one pair of sequences*/
			for(mpos=0,startpos[0]=0,startpos[1]=s_col[p]-1;mpos<msize;)
			{
				for (k=0;k<2;k++,mpos++)
				{
					marray[mpos].score=d1[startpos[0]][startpos[1]];
					marray[mpos].x=startpos[0];
					marray[mpos].y=startpos[1];
					for(j=startpos[0]+1,q=startpos[1]+1;j<s_col[i] && q<s_col[p];j++,q++)
					{
						if(d2[j][q] >= marray[mpos].score)/*we follow the principle of the righter the better*/
						{
							marray[mpos].score=d2[j][q];
							marray[mpos].x=j;
							marray[mpos].y=q;
						}
					}
					d2[marray[mpos].x][marray[mpos].y]=0;
				}
				/*shift the sequence*/
				if(startpos[1]>0) startpos[1]--;
				else 
				{		
					if (i==p)  break; /*avoid duplication when one sequence compare with itself*/
					else   startpos[0]++;
				}
			}
			/* select po->FIRST max segments from each pair of edges and refresh the peak matrix*/
			for(j=0,k=0;j<po->FIRST;j++,k=0)
			{
				maxindex=0;
				maxscore=0/*marray[maxindex].score*/;
				for(mpos=1;mpos<msize;mpos++)
				{
					if((maxscore < marray[mpos].score)) 
					{
						maxscore=marray[mpos].score;
						maxindex=mpos;
						k++;
					}
				}
				if (k==0) break;
				matrix[i][marray[maxindex].x]++;
				matrix[p][marray[maxindex].y]++;
				match[marray[maxindex].x]=TRUE;
				match1[marray[maxindex].y]=TRUE;
				marray[maxindex].score=0;
			}
		}
	}
	for(i=0;i<s_cols;i++)
        {
                free(d1[i]);
                free(d2[i]);
        }
        free(d1);
        free(d2);
}
/*************************************************************************************************************************/
static void pairwise_comparison_second ( bool **matrix3, int **matrix, int **matrix1, bool *match, bool *match1, int lower, int upper, struct dispos *marray, struct edgepos *marray1, int mpos, int msize, int maxscore, int maxindex, int *startpos)
{
	
	int i,j,p,q,k,j1=0,q1=0,min,min1;
	/*uglyTime("compare second start", i);*/
	continuous **d1, **d2;
        d1 = alloc2dd (s_cols,s_cols);
        d2 = alloc2dd (s_cols,s_cols);
	
	for(i=0;i<s_rows;i++)
	{
		/*uglyTime("second for i=%d", i);*/
		if ((i+1)%10==0)
		verboseDot();
		for(p=i;p<s_rows;p++)
		{
			/*initialize the values of d1 and d2*/
			for (j=0;j<s_col[i];j++)
				for (k=0;k<s_col[p];k++) 
			       		{d1[j][k]=0; d2[j][k]=0;}
			/*store values of d2 for the ith sequence of second alignment*/
			for(j=0;j<s_col[i];j++)
			{
				/*check matrix_no_continuous_equal*/
				if (matrix3[i][j])  continue;
				/*find local optimization of j in the first peak matrix*/
				min = MIN (j+po->local2, s_col[i]);
				j1=j; /*j1 is used to store the optimal position of current j in matrix*/
				min1 = matrix[i][j1];
				for (k=j-2;k<j+3;k++)
					if((k>=0&&k<=s_col[i])){
					if (matrix[i][k] >= min1) { j1=k; min1 = matrix[i][k];}
					}
				for(q=0;q<s_col[p];q++)
				{
					if (matrix3[p][q]) continue; /*check matrix_no_continuous_equal*/
					if (q==j && i==p)  {d1[j][q]=0; continue;}

					/*if (j>0 && q>0 && !matrix3[p][q-1] && !matrix3[i][j-1]) d1[j][q] = d1[j-1][q-1] - fre_matrix[p][seq_matrix[i][j-1]][q-1] + fre_matrix[p][seq_matrix[i][j+po->MOTIFLENGTH-1]][q+po->MOTIFLENGTH-1];
					else 
					{*/
						d1[j][q]=get_similarity_between_two_patterns (i,p,j,q,po->MOTIFLENGTH);
					/*}*/
					/*set a threshold for d1*/
					if (d1[j][q] ==0.0)  continue;
						d2[j][q] = improve_similarity_between_two_patterns(i, p, j, q, lower, upper, d1[j][q]);
					/*set a threshold for d2*/
					if (d2[j][q] == 0.0) { d2[j][q]=0; continue;}

					/*find local optimization of q in the first peak matrix*/
					min = MIN (q+po->local2, s_col[p]);
					q1=q; /*q1 is used to store the optimal position of q in matrix*/
					min1 = matrix[p][q1];
					for (k=q-2;k<q+3;k++)
						if((k>=0&&k<=s_col[p])){
						if (matrix[p][k] >= min1) { q1=k; min1 = matrix[p][k]+matrix[p][k+1];}
						}
					/* delete the edge with high weight but low value in matrix*/
					/*if (matrix[i][j1]+matrix[p][q1] <= 1) {d2[j][q]=0; continue;}*/
					/* calculate d2 base on d1*/
					d2[j][q]=d2[j][q]*(matrix[i][j1]+matrix[p][q1]);
				}	
			}
			/*select top2 in each alignment of one pair of sequences*/
			msize = 2*(s_col[i]+s_col[p]-1)-2;
			for(mpos=0,startpos[0]=0,startpos[1]=s_col[p]-1;mpos<msize;)
			{
				for(k=0;k<2;k++,mpos++)
				{
					marray[mpos].score=d2[startpos[0]][startpos[1]];
					marray[mpos].x=startpos[0];
					marray[mpos].y=startpos[1];
					for(j=startpos[0]+1,q=startpos[1]+1;j<s_col[i] && q<s_col[p];j++,q++)
						if(d2[j][q] >= marray[mpos].score)/*we follow the principle of the righter the better*/
						{
							marray[mpos].score=d2[j][q];
							marray[mpos].x=j;
							marray[mpos].y=q;
						}
					d2[marray[mpos].x][marray[mpos].y]=0;
				}	
				/*shift the sequence*/
				if(startpos[1]>0) startpos[1]--;
				else  
				{	
					if (i==p) break; 
					else startpos[0]++;
				}
			}
			/* select po->SECOND max segments and make the edges have no adjacent vertices, i.e. they are match*/
			for (j=0;j<s_cols;j++) { match[j]=FALSE; match1[j]=FALSE;}
			if (i==p)  q=msize/2+1;
			else  q=msize;
			for(j=0,k=0;j<po->SECOND;j++,k=0)
			{
				maxindex=0;
				maxscore=0;
				for(mpos=1;mpos<q;mpos++)
				{
					if((maxscore >= marray[mpos].score) || (match[marray[mpos].x]) || (match1[marray[mpos].y]) || (match[marray[mpos].x+1])||(match1[marray[mpos].y+1])) continue;
					if ((marray[mpos].x>0)&& (match[marray[mpos].x-1])) continue;
					if ((marray[mpos].y>0) && (match1[marray[mpos].y-1])) continue;
					maxscore=marray[mpos].score;
					maxindex=mpos;
					k++;
				}
				if (k==0) break;
				/*update matrix1 base on selected edges above*/
				matrix1[i][marray[maxindex].x]++;
				matrix1[p][marray[maxindex].y]++;
				match[marray[maxindex].x]=TRUE;
				match1[marray[maxindex].y]=TRUE;
				marray[maxindex].score=0;
			}
		}
	}
	for(i=0;i<s_cols;i++)
        {
                free(d1[i]);
                free(d2[i]);
        }
        free(d1);
        free(d2);
}
/***********************************************************************/
static void  pairwise_comparison_third (bool **matrix3, int **matrix, int **matrix1, bool *match, bool *match1, int lower, int upper, struct dispos *marray, struct edgepos *marray1, int mpos, int msize, int maxscore, int maxindex, int *startpos, int curelement1)
{
	
	/*int temp_length = po->MOTIFLENGTH;*/
	int row,col,min,min1,q1=0,j1=0, i,j,p,q,k;
	/*uglyTime("compare third start", i);*/
	/* get the 0-1 matrix base on the three times pairwise comparison, saved in acc_c; arr_c1 save the relative motif starting positions*/
	int ver = s_rows*po->TOPVERTICES;
	arr_c = alloc2d(ver,ver);
	arr_c1 = alloc2d(ver,ver);

	continuous **d1, **d2, **d3;
        d1 = alloc2dd (s_cols,s_cols);
        d2 = alloc2dd (s_cols,s_cols);
        d3 = alloc2dd (s_cols,s_cols);

	struct dispos **d4;
	d4=(struct dispos**)malloc(sizeof(struct dispos*)*s_cols);
        for(i=0;i<s_cols;i++)
                d4[i]=(struct dispos*)malloc(sizeof(struct dispos)*s_cols);

	/*initialize the VVmatrix arr_c and the accompany matrix arr_c1 which store the relative position of each pair of edges*/
	for (row = 0; row < ver; row++)
	{
		for (col = 0; col < ver; col++)
		{
			arr_c[row][col] = 0;
			arr_c1[row][col] = 0;
		}
		arr_c[row][row]=1;/*we set 1 to the elements on diagonal so that we can find clique base on arr_c*/
	}

	for (j=0;j<s_cols;j++)
		for (k=0;k<s_cols;k++) 
			{ d4[j][k].score=0; d3[j][k]=0;}

	for (i=0;i<s_rows;i++)
		for (j=0;j<s_cols;j++)
			matrix[i][j]=0;
	/*intiailze the matrix to store the peak value in the third alignment*/
	for (i=0;i<msize;i++)
		marray[i].score=0;
	for(i=0,curelement1=0;i<s_rows;i++)
	{
		/*uglyTime("third for i=%d", i);*/
		if ((i+1)%10==0) verboseDot();
		for(p=i;p<s_rows;p++)
		{
			/*initialize the values of d1, d4 and d3*/
			for (j=0;j<s_col[i];j++)
				for (k=0;k<s_col[p];k++)
					{ d4[j][k].score=0; d1[j][k]=0;	d3[j][k]=0; d2[j][k]=0;}
			for (j=0;j<s_cols;j++)
				{ match[j]=FALSE; match1[j]=FALSE; }

			for(j=0;j<s_col[i];j++)
			{
				if (matrix3[i][j])    continue; /*check matrix_no_continuous_equal*/
				/* find local optimization for j*/
				min = MIN (j+po->local3, s_col[i]);
				j1=j;
				min1 = matrix1[i][j1];
				for (k=j-2;k<j+3;k++)
					if((k>=0&&k<=s_col[i])){
					if (matrix1[i][k] > min1)
						{ j1=k; min1 = matrix1[i][k];}
					}
				for(q=0;q<s_col[p];q++)
				{
					if (matrix3[p][q])  continue;
					if (q==j && i==p) { d4[j][q].score = d1[j][q]=0; continue;}

					/*if (j>0 && q>0 && !matrix3[p][q-1] && !matrix3[i][j-1]) 
					{
						d1[j][q] = d1[j-1][q-1] - fre_matrix[p][seq_matrix[i][j-1]][q-1] + fre_matrix[p][seq_matrix[i][j+po->MOTIFLENGTH-1]][q+po->MOTIFLENGTH-1];
					}
					else 
					{*/
						d1[j][q]=get_similarity_between_two_patterns (i,p,j,q,po->MOTIFLENGTH);
						
					/*}*/
					if (d1[j][q] ==0.0)
					{
						if ((match[j]) && (match1[q])) continue;
						else { d4[j][q].score=0; continue;}
					}
					else
						d2[j][q] = improve_similarity_between_two_patterns(i, p, j, q, lower, upper, d1[j][q]);
					if (d2[j][q] ==0.0)
					{
						if ((match[j]) && (match1[q])) continue;
						else { /*d4[j][q].score=0;*/ continue;}
					}
					/*local optimize for d2*/
					min = MIN (q+po->local3, s_col[p]);
					q1=q;
					min1 = matrix1[p][q1];
					for (k=q-2;k<q+3;k++)
						if(k>=0&&k<=s_col[p]){
							if (matrix1[p][k] > min1)
								{ q1=k;	min1 = matrix1[p][k];}
						}

					if (matrix1[i][j1]+matrix1[p][q1] <=3)
					{
						/*d4[j][q].score =0;*/
						continue;
					}
					min = j1+q-j-q1;
					double log_height;
					log_height = d2[j][q]*log(matrix1[i][j1]+matrix1[p][q1]+1);
					if ((match[j1]) && (match1[q1]))
					{
						if ( (log_height > d4[j][q].score))
						{
							d4[j][q].score=log_height;
							d3[j1][q1]=min;
							d4[j][q].x=j1;
							d4[j][q].y=q1;
						}
						continue;
					}
					else
					{
						d4[j][q].x=j1;
						d4[j][q].y=q1;
						d4[j][q].score=log_height;
						d3[j1][q1]=min;
						match[j1]=TRUE;
						match1[q1]=TRUE;
					}
				}
			}
			/*select top2 in each alignment of one pair of sequences*/
			int k1=0,k2=0;
			msize = 2*(s_col[i]+s_col[p]-1)-2;
			for(mpos=0,startpos[0]=0,startpos[1]=s_col[p]-1;mpos<msize;)
			{
				for(k=0;k<2;k++,mpos++)
				{
					marray[mpos].score=d4[startpos[0]][startpos[1]].score;
					marray[mpos].x=startpos[0];
					marray[mpos].y=startpos[1];
					for(j=startpos[0],q=startpos[1];j<s_col[i] && q<s_col[p];j++,q++)
					{
						if(d4[j][q].score > marray[mpos].score )
						{
							marray[mpos].score=d4[j][q].score;
							marray[mpos].x=d4[j][q].x;
							marray[mpos].y=d4[j][q].y;
							marray[mpos].z=d3[d4[j][q].x][d4[j][q].y];
							k1=j;
							k2=q;
						}
					}
					d4[k1][k2].score=0;
				}
				if(startpos[1]>0) startpos[1]--;
				else
				{
					if (i==p)  break;
					else  startpos[0]++;
				}
			}
			/* select po->SECOND max segments*/
			bool *match_g;
			int local = floor (s_cols/3);
			AllocArray(match_g, s_cols);
			for (j=0; j<s_cols;j++)
				match_g[j] = TRUE;
			if (po->IS_global) 
			{
				for (j=local; j<(s_cols-local); j++)
					match_g[j] = FALSE;
			}
			else if (po->IS_local)
			{
				for (j=0; j<local; j++)
                                        match_g[j] = FALSE;
				for (j=s_cols-local; j<s_cols; j++)
                                        match_g[j] = FALSE;
			}
			else 
			{
				for (j=0; j<s_cols;j++)
	                                match_g[j] = FALSE;
			}
			for (j=0;j<s_cols;j++)
				{ match[j]=FALSE; match1[j]=FALSE;}

			int q1=0;
			if (i==p)  q1=po->THIRD/2;
			else  q1=po->THIRD;
			if (i==p) q=msize/2+1;
			else  q=msize;
			for(j=0,k=0;j<q1;j++,curelement1++,k=0)
			{
				maxindex=0;
				maxscore=0;
				marray1[curelement1].score=marray[maxindex].score;
				for(mpos=1;mpos<q;mpos++)
				{
					if ((match_g[marray[mpos].x])&&(match_g[marray[mpos].y])) continue;
					if((maxscore >= marray[mpos].score) || (match[marray[mpos].x]) || (match1[marray[mpos].y])) continue;
					maxscore=marray[mpos].score;
					maxindex=mpos;
					marray1[curelement1].score=marray[mpos].score;
					marray1[curelement1].x = marray[mpos].x;
					marray1[curelement1].y = marray[mpos].y;
					marray1[curelement1].x1 = i;
					marray1[curelement1].y1 = p;
					marray1[curelement1].z = marray[mpos].z;
					k++;
				}
				if (k==0) break;
				matrix[i][marray[maxindex].x]++;
				matrix[p][marray[maxindex].y]++;
				match[marray[maxindex].x]=TRUE;
				match1[marray[maxindex].y]=TRUE;
				marray[maxindex].score=0;
			}
		}
	}
	for(i=0;i<s_cols;i++)
        {
                free(d1[i]);
                free(d2[i]);
                free(d4[i]);
                free(d3[i]);
        }
        free(d1);
        free(d2);
        free(d4);
	free(d3);
}

/***********************************************************************/
static void select_topvertex (int **matrix, bool **matrix2, bool **matrix3)
/*select po->TOPVERTICES from each sequences so that we can construct a induced subgraph, saved in profile1*/
{
	int i,j,p,k, maxindex, mpos;
	for (i=0;i<s_rows;i++)
                for (j=0;j<s_cols;j++)
                        matrix2[i][j]=FALSE;
        int ver1=0;
        for (i=0;i<s_rows;i++)
        {
                /*consider sum weight  of two adjaent positions*/
                for (p=0;p<s_col[i];p++)
                        matrix[i][p]=matrix[i][p]+matrix[i][p+1];

                for (k=0;k<po->TOPVERTICES;k++)
                {
                        profile1[ver1].score=ver1;
                        maxindex=0; mpos=0; p=0;
                        if ((matrix[i][0])>0 && (!matrix2[i][1]) && (!matrix2[i][0]))
                                maxindex=matrix[i][0];
                        for (j=1;j<s_col[i];j++)
                        {
                                if ((matrix[i][j] < maxindex) || (matrix[i][j] < matrix[i][j-1]) || (matrix[i][j] < matrix[i][j+1]) || (matrix2[i][j]) || (matrix2[i][j-1]) || (matrix2[i][j+1]) || (matrix3[i][j]))  continue;
                                if (j<s_col[i]-1 && matrix[i][j] == matrix[i][j+1] && (matrix[i][j+2] >= matrix[i][j+1]||matrix[i][j-1] == matrix[i][j])) continue;
                                maxindex = matrix[i][j];
                                mpos=j;
                                p++;
                        }
                        if (p==0)
                                for (j=0;j<s_col[i];j++)
                                        if (!matrix2[i][j])
                                                { mpos=j; maxindex=matrix[i][j]; break;}
                        matrix2[i][mpos]=TRUE;
                        matrix[i][mpos]=-k-1;
                        profile1[ver1].x=i;
                        profile1[ver1].y=mpos;
                        ver1++;
                }
        }
}
/***********************************************************************/
static void get_final_graph (int msize1, int **matrix,  bool **matrix2,struct edgepos *marray1)
{
	int i,j,k1=0,k2=0,curelement1;
	/*update matrix2 so that the edges have not same vertex*/
        for (i=0;i<s_rows;i++)
                for (j=0;j<s_col[i];j++)
                        if (matrix2[i][j])
                                { matrix2[i][j+1]=TRUE; j++;}

        /*trace back for the toppo->THIRD edges so that we can get matrix arr_c and arr_c1*/
        for (curelement1=0;curelement1<msize1;curelement1++)
        {
                if ((matrix2[marray1[curelement1].x1][marray1[curelement1].x])&&(matrix2[marray1[curelement1].y1][marray1[curelement1].y]))
                {
                        if (marray1[curelement1].x==0)
                        {
                                i=-(matrix[marray1[curelement1].x1][marray1[curelement1].x])-1;
                                k1 = 1;
                        }
                        else
                        {
                                i=- MIN(matrix[marray1[curelement1].x1][marray1[curelement1].x],matrix[marray1[curelement1].x1][marray1[curelement1].x-1])-1;
                                if (matrix[marray1[curelement1].x1][marray1[curelement1].x]<matrix[marray1[curelement1].x1][marray1[curelement1].x-1])
                                        k1=1;
                                else
                                        k1=0;
                        }
                        if (marray1[curelement1].y==0)
                        {
                                j=-(matrix[marray1[curelement1].y1][marray1[curelement1].y])-1;
                                k2 = 1;
                        }
                        else
                        {
                                j=- MIN(matrix[marray1[curelement1].y1][marray1[curelement1].y],matrix[marray1[curelement1].y1][marray1[curelement1].y-1])-1;
                                if (matrix[marray1[curelement1].y1][marray1[curelement1].y]<matrix[marray1[curelement1].y1][marray1[curelement1].y-1])
                                        k2=1;
                                else
                                        k2=0;
                        }
                        arr_c[po->TOPVERTICES*marray1[curelement1].x1+i][po->TOPVERTICES*marray1[curelement1].y1+j]=1;
                        arr_c[po->TOPVERTICES*marray1[curelement1].y1+j][po->TOPVERTICES*marray1[curelement1].x1+i]=1;
                        if (k1==1 && k2==0)
                        {
                                arr_c1[po->TOPVERTICES*marray1[curelement1].x1+i][po->TOPVERTICES*marray1[curelement1].y1+j]= marray1[curelement1].z+1;
                                arr_c1[po->TOPVERTICES*marray1[curelement1].y1+j][po->TOPVERTICES*marray1[curelement1].x1+i]= marray1[curelement1].z+1;
                                continue;
                        }
                        if (k1==0 && k2==1)
                        {
                                arr_c1[po->TOPVERTICES*marray1[curelement1].x1+i][po->TOPVERTICES*marray1[curelement1].y1+j]= marray1[curelement1].z-1;
                                arr_c1[po->TOPVERTICES*marray1[curelement1].y1+j][po->TOPVERTICES*marray1[curelement1].x1+i]= marray1[curelement1].z-1;
                                continue;
                        }
                        arr_c1[po->TOPVERTICES*marray1[curelement1].x1+i][po->TOPVERTICES*marray1[curelement1].y1+j]= marray1[curelement1].z;
                        arr_c1[po->TOPVERTICES*marray1[curelement1].y1+j][po->TOPVERTICES*marray1[curelement1].x1+i]= marray1[curelement1].z;
                }
        }
}
/***********************************************************************/
/* the main signal strengthen part by compare sequences */
void compare_sequences(char **sequences)
{
	int i,p,j,curelement1=0;
	uglyTime("compare whole start", i);
	AllocArray (IsLengthEnough, s_rows);
	for (i=0; i<s_rows; i++)
		IsLengthEnough[i] = TRUE;
	/* set some advanced parameters for input motif length */
	parameter_choice();
	/* get the sequence length for each promoter which may be different*/
	for(i=0;i<s_rows;i++)
	{
		/*renew the length of s_col[] and s_cols when [L,U] and three patterns*/
		s_col[i]=strlen(sequences[i]);
		if (strlen(sequences[i])>s_cols)
                        s_cols=strlen(sequences[i]);

		if (s_col[i] < po->MOTIFLENGTH )
		{
		       	for (j=s_col[i];j<po->MOTIFLENGTH;j++)
			{
				sequences[i][j] = 'N';
			}
			s_col[i] = strlen(sequences[i]);
		}

		sequences[i][s_col[i]]='\0';
		s_col[i]=s_col[i]-(po->MOTIFLENGTH)+1;
	}
	/* we set s_cols larger than s_col[i] by 1 because we may need the memory in matrix2*/
	s_cols=s_cols-(po->MOTIFLENGTH)+2;
	/*build matrix*/
	int **matrix,**matrix1;/*the three peak matrix of the three alignment*/
	bool **matrix2,**matrix3;/*to select topvertices in the third peak matrix*/
	matrix=(int**)malloc(sizeof(int*)*s_rows);
	matrix1=(int**)malloc(sizeof(int*)*s_rows);
	matrix2=(bool**)malloc(sizeof(bool*)*s_rows);
	matrix3=(bool**)malloc(sizeof(bool*)*s_rows);
	for(i=0;i<s_rows;i++)
	{
		matrix[i]=(int*)malloc(sizeof(int)*s_cols);
		matrix1[i]=(int*)malloc(sizeof(int)*s_cols);
		matrix2[i]=(bool*)malloc(sizeof(bool)*s_cols);
		matrix3[i]=(bool*)malloc(sizeof(bool)*s_cols);
	}
	/* initialize the matrix and matrix1*/
	for(i=0;i<s_rows;i++)
		for(p=0;p<s_cols;p++) {	matrix[i][p]=0;	matrix1[i][p]=0;}

	/*initailize the matrix3 and update if there exist the continuous equal subsequence, e.g. AAAAA*/
	for (i=0;i<s_rows;i++)
		for (j=0;j<s_cols;j++)	matrix3[i][j]=FALSE;
	matrix3 = matrix_no_continuous_equal(matrix3);

	/*We need to select edges from each pair of sequences which compose a match graph base on match and match1*/
	bool *match, *match1;
	AllocArray(match, s_cols);
	AllocArray(match1,s_cols);
	
	/* define the two ends of motif */
	int upper,lower;
        double length;
        length=po->MOTIFLENGTH;
        lower = ceil(length/3);
        upper = po->MOTIFLENGTH-floor(length/3);

	/* generate matrix*/
        struct dispos *marray;
        struct edgepos *marray1;
        int msize,mpos=0,maxindex=0,maxscore=0,msize1,msize2;
	        /*select two maximum match for each location*/
        msize = 2*(2*s_cols-3)-2;
        msize1 = po->THIRD*s_rows*(s_rows+1)/2;/*store the topedges to generate the VVmatrix*/
        msize2 = s_rows*po->TOPVERTICES;/*store the information of topvertices in the third peak matrix*/
        marray=(struct dispos*)malloc(sizeof(struct dispos)*msize);
        marray1=(struct edgepos*)malloc(sizeof(struct edgepos)*msize1);
        profile1=(struct dispos*)malloc(sizeof(struct dispos)*msize2);

	/*initialize the marray1*/
	for (i=0; i<msize1; i++)
	{
		marray1[i].x=0;
		marray1[i].x1=0;
		marray1[i].y=0;
		marray1[i].y1=0;
	}
	for (i=0; i<msize; i++)
	{
		marray[i].x=0;
		marray[i].y=0;
	}

	int *startpos;/*store the two position of an edge*/	
	AllocArray (startpos,2);

	/* the first pairwise comparison between sequences*/
	pairwise_comparison_first (matrix3, matrix, match, match1, lower, upper,marray,marray1, mpos, msize, maxscore, maxindex, startpos);

	/* the second pairwise comparison between sequences*/
	pairwise_comparison_second (matrix3, matrix, matrix1, match, match1, lower, upper,marray,marray1, mpos, msize, maxscore, maxindex, startpos);

	/* the third passwise comparison between sequences*/
	pairwise_comparison_third (matrix3, matrix, matrix1, match, match1, lower, upper,marray,marray1, mpos, msize, maxscore, maxindex, startpos,curelement1);
	
	/*select po->TOPVERTICES from each sequences so that we can construct a induced subgraph, saved in profile1*/
	select_topvertex (matrix, matrix2, matrix3);	
	
	/*get the final graph on which we can find clique */	
	get_final_graph (msize1, matrix,  matrix2, marray1);
	
	printf ("\n");

	/*release the space*/
	if (po->no_enhance && po->MOTIFLENGTH == po->Up)
	{
		free(marray); free(marray1); free(match);free(match1);
		free(s_col);
		for(i=0;i<s_rows;i++)
		{
			free(matrix[i]); free(matrix1[i]); free(matrix2[i]);
		}
		free(matrix); free(matrix1); free(matrix2);
	}
	/*else
	{
		free(marray); free(marray1); free(match);free(match1);
		free(s_col);
		for(i=0;i<s_rows;i++)
		{
			free(matrix[i]); free(matrix1[i]); free(matrix2[i]);
		}
		free(matrix); free(matrix1); free(matrix2);
	}*/
	if (po->middle_enhance) uglyTime("Signal of middle part strengthening in %d sequences complete ", s_rows);
	else if (po->no_enhance) uglyTime("Signal of whole segment strengthening in %d sequences complete ", s_rows);
	else uglyTime("Signal of two ends strengthening in %d sequences complete ", s_rows);
}

/***********************************************************************/
