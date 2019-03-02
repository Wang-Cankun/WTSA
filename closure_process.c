/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 16, 2010
 * all the related process of closures are stored in the file        
 */


#include "closure_process.h"
static const int HEAP_SIZE = 3000000;
/************************************************************************/
discrete **frequency_matrix (char **sequence_temp, int a, int b, int motif_number)
/*get the frequency matrix used in caculating profile matrix*/
{
	int i,j=0,k;
	discrete **frequency;
        frequency = alloc2d (a,b);
        for (i=0; i<a; i++)
                for (j=0; j<b; j++)
                        frequency[i][j]=0;
/*printf ("%d\t%d\t%d\n",a,b,motif_number);*/
        for (i=0; i<motif_number; i++)
                for (k=0;k<b;k++)
                {
                        if (sequence_temp[i][k]=='A' || sequence_temp[i][k]=='a') frequency[1][k]=frequency[1][k]+1;
                        else if (sequence_temp[i][k]=='G' || sequence_temp[i][k]=='g') frequency[2][k]=frequency[2][k]+1;
                        else if (sequence_temp[i][k]=='C' || sequence_temp[i][k]=='c') frequency[3][k]=frequency[3][k]+1;
                        else if (sequence_temp[i][k]=='T' || sequence_temp[i][k]=='t') frequency[4][k]=frequency[4][k]+1;
                }
	return frequency;
}
/************************************************************************/
continuous **get_profile (char **sequence_temp, int a, int b, int motif_number)
/*get the profile matrix of a closure*/
{
	continuous **scoreM;
	discrete **frequency;
	int i,j;
	scoreM = alloc2dd (a,b);
/*	scoreM = (continuous**)malloc(sizeof(continuous*)*a);
        for(i=0;i<a;i++)
                scoreM[i]=(continuous*)malloc(sizeof(continuous)*b);*/

	for (i=0; i<a; i++)
		for (j=0; j<b; j++)
			scoreM[i][j]=0;
	frequency = frequency_matrix (sequence_temp, a, b, motif_number);


	for (i=1;i<a;i++)
		for (j=0;j<b;j++)
			scoreM[i][j]= log2((frequency[i][j]+p_markov[i][0])/(p_markov[i][0]*(motif_number+1)));
	for (i=0;i<a;i++)
		free (frequency[i]);
	free (frequency);
	return scoreM;
}

/************************************************************************/
continuous get_similarity_between_two_patterns(int seq1, int seq2, int pos1, int pos2, int motif_length)
{
        int i=0;
        continuous num=0;
	/*continuous binomial18[] ={
0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,2.2
,2.9
,3.6
,4.4
,5.4
,6.4
,7.6
,9.1
,10.8};*/

continuous binomial14[] ={
0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,2.7
,3.5
,4.4
,5.5
,6.8
,8.4};
/*
 for(i=0;i<motif_length;++i){
	if(pos1+i<MAX_SEQUENCE_LENGTH){
		num+=fre_matrix[seq2][seq_matrix[seq1][pos1+i]][pos2+i];
	}
}
*/
num = fre_matrix[seq2][seq_matrix[seq1][pos1+1]][pos2+1] 
+ fre_matrix[seq2][seq_matrix[seq1][pos1+2]][pos2+2]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+3]][pos2+3]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+4]][pos2+4]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+5]][pos2+5]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+6]][pos2+6]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+7]][pos2+7]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+8]][pos2+8]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+9]][pos2+9]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+10]][pos2+10]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+11]][pos2+11]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+12]][pos2+12]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+13]][pos2+13]
+ fre_matrix[seq2][seq_matrix[seq1][pos1+14]][pos2+14];

if(num < 8){
	return 0.0;
	
}else {
        return (binomial14[(int)num]);
	}
}

/************************************************************************/
continuous improve_similarity_between_two_patterns(int seq1, int seq2, int pos1, int pos2, int lower, int upper, continuous init)
{
       int pos=0;
      /* double binomial18[] ={
0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,2.2
,2.9
,3.6
,4.4
,5.4
,6.4
,7.6
,9.1
,10.8};*/

double binomial14[] ={
0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,0.0
,2.7
,3.5
,4.4
,5.5
,6.8
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4
,8.4};
       continuous num= init;
       if (po->middle_enhance)
       {
               for (pos=lower;pos<upper;pos++)
		       num += fre_matrix[seq2][seq_matrix[seq1][pos1+pos]][pos2+pos]*(po->end_weight-1); 
       }
       if (po->no_enhance)
               num = 2*num;
       else
       {
               /* we consider the conserve proporty of two ends of motif*/
               num=(po->end_weight)*num;
               for (pos=lower;pos<upper;pos++){
		       num -= fre_matrix[seq2][seq_matrix[seq1][pos1+pos]][pos2+pos]*(po->end_weight-1);
	       }
		       
       }
       /*if (po->SequenceWeight){
               num = num*(SequenceWeight[seq1]+SequenceWeight[seq2]);
       }*/

	
	/*printf("%d\t%d\t\n",max1,max2);*/
	/*num = num * (log10f((float)height_matrix[seq1][pos1])/log10f((float)max1)+log10f((float)height_matrix[seq2][pos2])/log10f((float)max2));*/
	/*printf("%d\t\n",height_matrix[seq1][pos1]);*/
	/*if(num<15.0){
       		num = (binomial14[(int)num]);
	}else {
	num = (binomial14[13]);
	}*/
	
	/*printf("%f\t%d\t\n",num,max2);*/
	/*num = num * (log1pf((float)height_matrix[seq1][pos1])/logf((float)max1)+log1pf((float)height_matrix[seq2][pos2])/logf((float)max2));
	
	if(max2==1){num=0.0;}*/
	
	/*printf("%f\t%d\t\n",num,max2);*/
	if(num < 8){
	return 0;
	}else {
        return (binomial14[(int)num]);
	}
}

/* C function to find maximum in arr[] of size n*/
int largest(int arr[], int n)
{
    int i;
    
    /* Initialize maximum element*/
    int max = arr[0];
 
    /* Traverse array elements from second and
     compare every element with current max  */
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
 
    return max;
}
 

/************************************************************************/
continuous palindromic_pro (char *consensus, int forepart)
/* get the proportion of palindromic of forepart */
{
	int i, length=0;
       	continuous palindromic_num=0, palindromic_max=0;
	continuous proportion;
	length = strlen(consensus);
	char *palindromic;
	AllocArray (palindromic, length);
	for (i=0;i<length;i++)
	{
		if (consensus[i]=='A'|| consensus[i]=='a') { palindromic[length-i-1]='T'; consensus[i]='A';}
		if ((consensus[i]=='T')||(consensus[i]=='t')) {palindromic[length-i-1]='A';consensus[i]='T';}
		if ((consensus[i]=='C')||(consensus[i]=='c')) {palindromic[length-i-1]='G';consensus[i]='C';}
		if ((consensus[i]=='G')||(consensus[i]=='g')) {palindromic[length-i-1]='C';consensus[i]='G';}
	}
	/*printf ("1: %d\t%s\n",length,palindromic);*/
	for (i=0;i<forepart;i++)
		if (consensus[i] == palindromic[i]) palindromic_num++;
	if (palindromic_num > palindromic_max)
	{
		palindromic_max = palindromic_num;
		palindromic_num = 0;
	}
	for (i=0;i<forepart-1;i++)
                if (consensus[i] == palindromic[i+1]) palindromic_num++;
	if (palindromic_num > palindromic_max)
        {
                palindromic_max = palindromic_num;
                palindromic_num = 0;
        }
        palindromic_num = 0;
	for (i=0;i<forepart;i++)
                if (consensus[i+1] == palindromic[i]) palindromic_num++;
	if (palindromic_num > palindromic_max) palindromic_max = palindromic_num;
	/*printf ("2: %s\n",palindromic);
	printf ("3: %f\t%d\n",palindromic_max,forepart);*/
	proportion = palindromic_max/forepart;
	return proportion;
}
/************************************************************************/
bool palindromic_pro_1 (char *consensus, int forepart)
/* get the proportion of palindromic of forepart */
{
	int i, length=0;
	int *num;
	AllocArray (num, 2);
	num[0] = num[1] =0;
       	continuous palindromic_num_7=0, palindromic_num_5=0;
	length = strlen(consensus);
	char *palindromic;
	AllocArray (palindromic, length);
	
	/* check left AA/TT and right GC in to-be-identify palindromic*/
	if (!(((consensus[po->Palin1]=='A'||consensus[po->Palin1]=='a')&&(consensus[po->Palin1+1]=='A'||consensus[po->Palin1+1]=='a')&&(consensus[length-po->Palin1-1]=='C'||consensus[length-po->Palin1-1]=='c')&&(consensus[length-po->Palin1-2]=='G'||consensus[length-po->Palin1-2]=='g')) || 
	      ((consensus[po->Palin1]=='G'||consensus[po->Palin1]=='g')&&(consensus[po->Palin1+1]=='C'||consensus[po->Palin1+1]=='c')&&(consensus[length-po->Palin1-1]=='T'||consensus[length-po->Palin1-1]=='t')&&(consensus[length-po->Palin1-2]=='T'||consensus[length-po->Palin1-2]=='t')))) 
	/* check left AA/TT and right CG in to-be-identify palindromic*/
	/*if (!(((consensus[po->Palin1]=='A'||consensus[po->Palin1]=='a')&&(consensus[po->Palin1+1]=='A'||consensus[po->Palin1+1]=='a')&&(consensus[length-po->Palin1-1]=='G'||consensus[length-po->Palin1-1]=='g')&&(consensus[length-po->Palin1-2]=='C'||consensus[length-po->Palin1-2]=='c')) || 
	      ((consensus[po->Palin1]=='C'||consensus[po->Palin1]=='c')&&(consensus[po->Palin1+1]=='G'||consensus[po->Palin1+1]=='g')&&(consensus[length-po->Palin1-1]=='T'||consensus[length-po->Palin1-1]=='t')&&(consensus[length-po->Palin1-2]=='T'||consensus[length-po->Palin1-2]=='t')))) */
	{
		return 0;
	}	

	/* construct palindromic sequence*/
	for (i=0;i<length;i++)
	{
		if (consensus[i]=='A'|| consensus[i]=='a') { palindromic[length-i-1]='T'; consensus[i]='A';}
		if ((consensus[i]=='T')||(consensus[i]=='t')) {palindromic[length-i-1]='A';consensus[i]='T';}
		if ((consensus[i]=='C')||(consensus[i]=='c')) {palindromic[length-i-1]='G';consensus[i]='C';}
		if ((consensus[i]=='G')||(consensus[i]=='g')) {palindromic[length-i-1]='C';consensus[i]='G';}
	}

	/*compare for the 7-length part*/
	for (i=0;i<po->Palin1;i++)
		if (consensus[i] == palindromic[i]) 
			palindromic_num_7++;
	if (palindromic_num_7 < po->Palin1-po->Palin3)
		return 0;

	/*compare for the 5-length part*/
	for (i=po->Palin1+2;i<po->Palin1+po->Palin2+2;i++)
		if (consensus[i] == palindromic[i]) 
			palindromic_num_5++;
	if(palindromic_num_5 < po->Palin2-po->Palin3)
		return 0;

	/*printf ("%.0f\t%.0f\n",palindromic_num_7,palindromic_num_5);*/	
	return 1;
}
/************************************************************************/
bool check_palindromic (char **sequence_palin, int a, int b, int motif_number,int forepart)
{	
	int i, max=0;
	discrete **frequency;
	continuous proportion=0;
	frequency = frequency_matrix (sequence_palin, a, b, motif_number);
	char *consensus;
	AllocArray(consensus,b);
	for (i=0;i<b;i++)
	{
		max=0;
		if(frequency[1][i] > max) {consensus[i]='A'; max = frequency[1][i];}
		if(frequency[2][i] > max) {consensus[i]='G'; max = frequency[2][i];}
		if(frequency[3][i] > max) {consensus[i]='C'; max = frequency[3][i];}	
		if(frequency[4][i] > max) {consensus[i]='T'; max = frequency[4][i];}
	}		
	proportion = palindromic_pro (consensus, forepart);
/*	printf ("%3.2f\n",proportion);*/
	if (proportion > 0.6) { return TRUE;}
	else {return FALSE;}
}
/************************************************************************/
continuous **get_palindromic_profile (continuous **scoreM, int length_local)
/*get the palindromic profile*/
{
        int j;
        continuous **scoreM_panlindromic;
        scoreM_panlindromic = alloc2dd (5,length_local);
        for (j=0;j<length_local;j++)
        {
                scoreM_panlindromic[0][length_local-j-1] = scoreM[0][j];
                scoreM_panlindromic[1][length_local-j-1] = scoreM[4][j];
                scoreM_panlindromic[4][length_local-j-1] = scoreM[1][j];
                scoreM_panlindromic[2][length_local-j-1] = scoreM[3][j];
                scoreM_panlindromic[3][length_local-j-1] = scoreM[2][j];
        }
        return scoreM_panlindromic;
}
/************************************************************************/

continuous mirror_pro (char *consensus, int forepart)
/* get the proportion of mirror of forepart */
{
	int i, length=0;
       	continuous mirror_num=0, mirror_max=0;
	continuous proportion;
	length = strlen(consensus);
	char *mirror;
	AllocArray (mirror, length);
/*	printf ("%s\t",consensus);*/
	for (i=0;i<length;i++)
	{
		if (consensus[i]=='A'|| consensus[i]=='a') { mirror[length-i-1]='A'; consensus[i]='A';}
		if ((consensus[i]=='T')||(consensus[i]=='t')) { mirror[length-i-1]='T';consensus[i]='T';}
		if ((consensus[i]=='C')||(consensus[i]=='c')) { mirror[length-i-1]='C';consensus[i]='C';}
		if ((consensus[i]=='G')||(consensus[i]=='g')) { mirror[length-i-1]='G';consensus[i]='G';}
	}
/*	printf ("%s\t",mirror);*/
	for (i=0;i<forepart;i++)
		if (consensus[i] == mirror[i]) mirror_num++;
	if (mirror_num > mirror_max)
	{
		mirror_max = mirror_num;
		mirror_num = 0;
	}
	for (i=0;i<forepart;i++)
                if (consensus[i] == mirror[i+1]) mirror_num++;
	if (mirror_num > mirror_max)
        {
                mirror_max = mirror_num;
                mirror_num = 0;
        }
	for (i=0;i<forepart;i++)
                if (consensus[i+1] == mirror[i]) mirror_num++;
	if (mirror_num > mirror_max) mirror_max = mirror_num;
/*	printf ("%f\t%d\n",mirror_max,forepart);*/
	proportion = mirror_max/forepart;
	return proportion;
}
/************************************************************************/
bool check_mirror (char **sequence_mirror, int a, int b, int motif_number,int forepart)
{	
	int i, max=0;
	discrete **frequency;
	continuous proportion=0;
	frequency = frequency_matrix (sequence_mirror, a, b, motif_number);
	char *consensus;
	AllocArray(consensus,b);
	for (i=0;i<b;i++)
	{
		max = 0;
		if(frequency[1][i] > max) {consensus[i]='A'; max = frequency[1][i];}
		if(frequency[2][i] > max) {consensus[i]='G'; max = frequency[2][i];}
		if(frequency[3][i] > max) {consensus[i]='C'; max = frequency[3][i];}	
		if(frequency[4][i] > max) {consensus[i]='T'; max = frequency[4][i];}
	}		
	proportion = mirror_pro (consensus, forepart);
/*	printf ("%3.2f\n",proportion);*/
	if (proportion > 0.6) { return TRUE;}
	else {return FALSE;}
}

/************************************************************************/
continuous **impovre_profle_mirror (continuous **scoreM, double  length_local_1,  char **sequence_1, int a ,int motif_number)
/*  improve the score matrix base on the conserve property of two ends of motif, considering mirror*/
{
	int forepart=floor(length_local_1/3),endpart=ceil(length_local_1/3),delIndex=0;
        continuous sum_scoreM=0;
        int i,j;
        for (j=0;j<forepart;j++)
        {
                scoreM[0][j]=0;
                for (i=1;i<5;i++)
                        if (scoreM[i][j] > scoreM[0][j]) scoreM[0][j]=scoreM[i][j];
                sum_scoreM = sum_scoreM + scoreM[0][j];
        }
        for (j=0;j<forepart;j++)
                if (scoreM[0][j]<((sum_scoreM*0.6)/forepart))
                {
                        if (delIndex==1)
                                for (i=1;i<5;i++)
                                        scoreM[i][j]=0;
                        else
                                delIndex=1;
                }
        sum_scoreM =0; delIndex=0;
        for (j=length_local_1-1;j>length_local_1-endpart-1;j--)
        {
                scoreM[0][j]=0;
                for (i=1;i<5;i++)
                        if (scoreM[i][j] > scoreM[0][j]) scoreM[0][j]=scoreM[i][j];
                sum_scoreM = sum_scoreM + scoreM[0][j];
        }
	for (j=length_local_1-1;j>length_local_1-endpart-1;j--)
                if (scoreM[0][j]<((sum_scoreM*0.6)/forepart))
                {
                        if (delIndex==1)
                                for (i=1;i<5;i++)
                                        scoreM[i][j]=0;
                        else
                                delIndex=1;
                }

	if (check_mirror(sequence_1, a, length_local_1,motif_number, forepart))
	{
		for (j=forepart;j<length_local_1-endpart;j++)
                	for (i=1;i<5;i++)
                        	scoreM[i][j]=0;
	}
	else
	{
        	for (j=forepart;j<length_local_1-endpart;j++)
                	for (i=1;i<5;i++)
				scoreM[i][j]=scoreM[i][j]/2; 
	}
        return scoreM;
}
/********************************************************************************************/
continuous **impovre_profle_palindromic (continuous **scoreM, double  length_local_1,  char **sequence_1, int a ,int motif_number)
/*  improve the score matrix base on the conserve property of two ends of motif, considering palindromic*/
{
	int forepart=floor(length_local_1/3),endpart=ceil(length_local_1/3),delIndex=0;
        continuous sum_scoreM=0;
        int i,j;
        for (j=0;j<forepart;j++)
        {
                scoreM[0][j]=0;
                for (i=1;i<5;i++)
                        if (scoreM[i][j] > scoreM[0][j]) scoreM[0][j]=scoreM[i][j];
                sum_scoreM = sum_scoreM + scoreM[0][j];
        }
        for (j=0;j<forepart;j++)
                if (scoreM[0][j]<((sum_scoreM*0.6)/forepart))
                {
                        if (delIndex==1)
                                for (i=1;i<5;i++)
                                        scoreM[i][j]=0;
                        else
                                delIndex=1;
                }
        sum_scoreM =0; delIndex=0;
        for (j=length_local_1-1;j>length_local_1-endpart-1;j--)
        {
                scoreM[0][j]=0;
                for (i=1;i<5;i++)
                        if (scoreM[i][j] > scoreM[0][j]) scoreM[0][j]=scoreM[i][j];
                sum_scoreM = sum_scoreM + scoreM[0][j];
        }
	for (j=length_local_1-1;j>length_local_1-endpart-1;j--)
                if (scoreM[0][j]<((sum_scoreM*0.6)/forepart))
                {
                        if (delIndex==1)
                                for (i=1;i<5;i++)
                                        scoreM[i][j]=0;
                        else
                                delIndex=1;
                }

	if (check_palindromic (sequence_1, a, length_local_1,motif_number, forepart))
	{
		for (j=forepart;j<length_local_1-endpart;j++)
                	for (i=1;i<5;i++)
                        	scoreM[i][j]=0;
	}
	else
	{
        	for (j=forepart;j<length_local_1-endpart;j++)
                	for (i=1;i<5;i++)
				scoreM[i][j]=scoreM[i][j]/2; 
	}
        return scoreM;
}
/************************************************************************/
continuous **impovre_profle (continuous **scoreM, double  length_local_1)
/*improve the score matrix base on the conserve property of two ends of motif*/
{
	int forepart=floor(length_local_1/3),endpart=ceil(length_local_1/3),delIndex=0;
        continuous sum_scoreM=0;
	int i,j;
        for (j=0;j<forepart;j++)
        {
                scoreM[0][j]=0;
                for (i=1;i<5;i++)
			if (scoreM[i][j] > scoreM[0][j]) scoreM[0][j]=scoreM[i][j];
                sum_scoreM = sum_scoreM + scoreM[0][j];
        }
        for (j=0;j<forepart;j++)
	        if (scoreM[0][j]<((sum_scoreM*0.6)/forepart))
                {
        	        if (delIndex==1)
                	        for (i=1;i<5;i++)
                        	        scoreM[i][j]=0;
			else
				delIndex=1;
                }
	sum_scoreM =0; delIndex=0;
        for (j=length_local_1-1;j>length_local_1-endpart-1;j--)
        {
        	scoreM[0][j]=0;
                for (i=1;i<5;i++)
                	if (scoreM[i][j] > scoreM[0][j]) scoreM[0][j]=scoreM[i][j];
                sum_scoreM = sum_scoreM + scoreM[0][j];
        }
	for (j=length_local_1-1;j>length_local_1-endpart-1;j--)
        	if (scoreM[0][j]<((sum_scoreM*0.6)/forepart))
                {
                	if (delIndex==1)
                        	for (i=1;i<5;i++)
                                	scoreM[i][j]=0;
			else
	                        delIndex=1;
		}
	if (po->end_weight > 1)
	{
		for (j=forepart;j<length_local_1-endpart;j++)
        		for (i=1;i<5;i++)
				scoreM[i][j]=scoreM[i][j]/2;  /*if plalindome, scoreM[i][j]=0*/
	}
	return scoreM;
}

/************************************************************************/
continuous get_IC (char **sequence_temp, int a, int b, int motif_number)
/*get the information content of a closure base on profile matrix*/
{
	int i,j;
	continuous IC = 0;
        continuous **scoreM;
        discrete **frequency;
	frequency = frequency_matrix (sequence_temp, a, b, motif_number);
	scoreM = get_profile (sequence_temp, a, b, motif_number);
        for (i=1;i<a;i++)
                for (j=0;j<b;j++)
			IC = IC + ((frequency[i][j]+p_markov[i][0])*scoreM[i][j])/(motif_number+1);
	for (i=0;i<a;i++)
        {
                free(frequency[i]);
                free(scoreM[i]);
        }
        free(frequency);
        free(scoreM);
        return IC;
}
/************************************************************************/
continuous IC_closures_1 (char **sequence_1, char **sequence_2, int a, int b, int motif_number_1, int motif_number_2, int p ,int q)
/*get information content of the first closure between two closures defined in our paper*/
{
	int i,j;
	int b_1 = b/2;
	continuous IC = 0;
	continuous **scoreM_1, **scoreM_2;
        discrete **frequency_1, **frequency_2;
	frequency_1 = frequency_matrix (sequence_1, a, b, motif_number_1);
	frequency_2 = frequency_matrix (sequence_2, a, b, motif_number_2);
	scoreM_1 = get_profile (sequence_1, a, b, motif_number_1);	
	scoreM_2 = get_profile (sequence_2, a, b, motif_number_2);
	for (i=1;i<a;i++)
                for (j=0;j<b_1;j++)
                        IC = IC + ((frequency_2[i][q+j]+p_markov[i][0])*scoreM_1[i][p+j])/(motif_number_2+1);
	for (i=0;i<a;i++)
        {
                free(frequency_1[i]);
                free(frequency_2[i]);
                free(scoreM_1[i]);
                free(scoreM_2[i]);
        }
        free(frequency_1);
        free(frequency_2);
        free(scoreM_1);
        free(scoreM_2);
        return IC;
}
/************************************************************************/
continuous IC_closures_2 (char **sequence_1, char **sequence_2, int a, int b, int motif_number_1, int motif_number_2, int p ,int q)
/*get information content of the second closure between two closures defined in our paper*/
{
	int i,j;
	int b_1 = b/2;
        continuous IC = 0;
        continuous **scoreM_1, **scoreM_2;
        discrete **frequency_1, **frequency_2;
        frequency_1 = frequency_matrix (sequence_1, a, b, motif_number_1);   
        frequency_2 = frequency_matrix (sequence_2, a, b, motif_number_2);
        scoreM_1 = get_profile (sequence_1, a, b, motif_number_1);   
        scoreM_2 = get_profile (sequence_2, a, b, motif_number_2);
        for (i=1;i<a;i++)
                for (j=0;j<b_1;j++)
                        IC = IC + ((frequency_1[i][p+j]+p_markov[i][0])*scoreM_2[i][q+j])/(motif_number_1+1);
	for (i=0;i<a;i++)
	{	
		free(frequency_1[i]);
		free(frequency_2[i]);
		free(scoreM_1[i]);
		free(scoreM_2[i]);
	}
	free(frequency_1);
        free(frequency_2);
        free(scoreM_1);
        free(scoreM_2);
        return IC;
}
/************************************************************************/
char **get_sequences_from_closures (Closures *aa)
{
	int i,j,k,genome_num=0;
	char **sequence_aa;
	bool flag = FALSE;
	for (i=0;i<reference_oper_num;i++)
	{
		if (sameString(genome[i]->oper_name,aa->name))
		{
			flag = TRUE;
			genome_num = i;
			break;
		}	
	}
        sequence_aa = alloc2c (aa->closure_rows,aa->length+extend_len);
	for (i=0; i<aa->closure_rows; i++)
	{
		k=0;
		for (j=dsItem(aa->position,i);j< (dsItem(aa->position,i)+aa->length+extend_len); j++)
		{
			if (flag) sequence_aa[i][k] = genome[genome_num]->sequences_r[dsItem(aa->sequence,i)][j];
			else sequence_aa[i][k] = sequences[dsItem(aa->sequence,i)][j];
			k++;
		}
	}
	return sequence_aa;
}

/************************************************************************/
char **get_2L_sequeces_from_closures (Closures *aa)
/*get sequences with 2L base on a closure*/
{
	int i,j,k, left, right, left_1, right_1, length1;
	char **sequence_aa;
	int genome_num = 0;
	bool flag = FALSE;
	for (i=0;i<reference_oper_num;i++)
        {
                if (sameString(genome[i]->oper_name,aa->name))
                {
                        flag = TRUE;
                        genome_num = i;
                        break;
                }
        }
	sequence_aa = alloc2c (aa->closure_rows,2*(aa->length+extend_len));
	for (i=0; i<aa->closure_rows; i++)
        {
                k=0;
                left_1 = dsItem(aa->position,i) - floor((aa->length+extend_len)/2);
                left = MAX(0,left_1);
		if (flag) length1= strlen(genome[genome_num]->sequences_r[dsItem(aa->sequence,i)])-1;
		else length1=strlen(sequences[dsItem(aa->sequence,i)])-1;
                right_1 = dsItem(aa->position,i) + aa->length+extend_len + ceil((aa->length+extend_len)/2);
                right = MIN(length1,right_1);
                if (left_1<0)
                {
                        for (j=0; j< ABS(left_1); j++)
                        {
                                sequence_aa[i][k] = 'N';
                                k++;
                        }
                }
                for (j=left; j< right; j++)
                {
                        if (flag) sequence_aa[i][k] = genome[genome_num]->sequences_r[dsItem(aa->sequence,i)][j];
			else sequence_aa[i][k] = sequences[dsItem(aa->sequence,i)][j];
                        k++;
                }
                if (right_1 > length1)
                        for (j=0;j<(right_1 - length1);j++)
                        {
                                sequence_aa[i][k] = 'N';
                                k++;
                        }
        }
	return sequence_aa;
}
/************************************************************************/
continuous similarity_closures(Closures *aa, Closures *bb)
/*get the similarity between two different closures*/
{
	int i,j;
	continuous IC_a, IC_b, IC_ab, IC_ba, IC_simi = 0, IC_temp =0;
	char **sequence_aa_1, **sequence_aa_2, **sequence_bb_1, **sequence_bb_2;
	/*do not re-alloc2c equence_aa_1 if it get from get_sequences_from_closures, otherwise memory will overflow*/
	sequence_aa_1 = get_sequences_from_closures (aa);
	sequence_aa_2 = get_2L_sequeces_from_closures (aa);
	sequence_bb_1 = get_sequences_from_closures (bb);
        sequence_bb_2 = get_2L_sequeces_from_closures (bb);
	IC_a = get_IC (sequence_aa_1, 5, aa->length+extend_len, aa->closure_rows);
	IC_b = get_IC (sequence_bb_1, 5, bb->length+extend_len, bb->closure_rows);
	for (i=0;i<aa->length+extend_len;i++)
	{
		for (j=0; j< bb->length+extend_len; j++)
		{
			IC_ab = IC_closures_1 (sequence_aa_2, sequence_bb_2, 5, 2*(aa->length+extend_len), aa->closure_rows, bb->closure_rows, i ,j);
			IC_ba = IC_closures_2 (sequence_aa_2, sequence_bb_2, 5, 2*(bb->length+extend_len), aa->closure_rows, bb->closure_rows, i ,j);
			IC_temp = (IC_ab+IC_ba)/(IC_a+IC_b);
			if (IC_temp > IC_simi)	IC_simi = IC_temp;
		}
	}
	for (i=0;i<aa->closure_rows;i++)
		free(sequence_aa_1[i]);
	free (sequence_aa_1);
	for (i=0;i<aa->closure_rows;i++)
                free(sequence_aa_2[i]);
        free (sequence_aa_2);
	for (i=0;i<bb->closure_rows;i++)
                free(sequence_bb_1[i]);
        free (sequence_bb_1);
	for (i=0;i<bb->closure_rows;i++)
                free(sequence_bb_2[i]);
        free (sequence_bb_2);
	return IC_simi;
}
/************************************************************************/
continuous *get_column_IC (Closures **aa, int clos)
{
        continuous *I_11;
        int p , q;
        int length = 2*(aa[clos]->length+extend_len);
        AllocArray (I_11, length);
        for (p = 0; p < length; p++)
        {
                I_11[p] = 0;
                for (q=1; q<5; q++)
                {
                        I_11[p] = I_11[p] + ((fre2_all[clos][q][p]+p_markov[q][0])*sco2_all[clos][q][p])/(aa[clos]->closure_rows+1);
                }
        }
        return I_11;
}
/************************************************************************/
continuous *get_column_IC_12 (Closures **aa, int clos_1, int clos_2, int startpos_1, int startpos_2)
{
        continuous *I_12;
        int p , q;
        int length = aa[clos_1]->length+aa[clos_2]->length+extend_len+extend_len;
        AllocArray (I_12, length);
        int pos_max = MAX (startpos_1,startpos_2);
        int length_new = length-pos_max;
        for (p = 0; p < length_new; p++)
        {
                I_12[p] = 0;
                for (q=1; q<5; q++)
                {
		        I_12[p] = I_12[p] + ((fre2_all[clos_1][q][startpos_1+p]+p_markov[q][0])*sco2_all[clos_2][q][startpos_2+p])/(aa[clos_1]->closure_rows+1);
        	}
		fflush(stdout);
	}
	fflush(stdout);
        return I_12;
}
/************************************************************************/
continuous *get_column_IC_21 (Closures **aa, int clos_1, int clos_2, int startpos_1, int startpos_2)
{
        continuous *I_21;
        int p , q;
        int length = aa[clos_1]->length+aa[clos_2]->length+extend_len+extend_len;
        AllocArray (I_21, length);
        int pos_max = MAX (startpos_1,startpos_2);
        int length_new = length-pos_max;
        for (p = 0; p < length_new; p++)
        {
                I_21[p] = 0;
                for (q=1; q<5; q++)
		{
                        I_21[p] = I_21[p] + ((fre2_all[clos_2][q][startpos_2+p]+p_markov[q][0])*sco2_all[clos_1][q][startpos_1+p])/(aa[clos_2]->closure_rows+1);
		}
		startpos_1=startpos_1; startpos_2=startpos_2;
		fflush(stdout);
        }
	fflush(stdout);
        return I_21;
}
/************************************************************************/
bool *clean_up_closures (Closures **aa, int closure_id, continuous threshold)
{
	bool *IS_duplicate; /*remove duplicate from the closures*/
        AllocArray (IS_duplicate,closure_id);
        int i,j,ii,jj,p,q;
        char **sequence_aa_1, **sequence_aa_2;
        
	for (i=0;i<closure_id;i++)
                IS_duplicate[i] = FALSE;
	
        continuous IC_1=0,IC_2=0,IC_12=0,IC_21=0,IC_temp =0,similarity=-1;
        fre_all = (discrete***)malloc(sizeof(discrete**)*closure_id);
        fre2_all = (discrete***)malloc(sizeof(discrete**)*closure_id);
        sco_all = (continuous***)malloc(sizeof(continuous**)*closure_id);
        sco2_all = (continuous***)malloc(sizeof(continuous**)*closure_id);
        for (i=0;i<closure_id;i++)
        {
                fre_all[i] = alloc2d (5, aa[i]->length+extend_len);
                fre2_all[i] = alloc2d (5, 2*aa[i]->length+extend_len);
                sco_all[i] = alloc2dd (5, aa[i]->length+extend_len);
                sco2_all[i] = alloc2dd (5, 2*aa[i]->length+extend_len);
        }
	for (i=0;i<closure_id;i++)
        {
		sequence_aa_1 = get_sequences_from_closures (aa[i]);
                sequence_aa_2 = get_2L_sequeces_from_closures (aa[i]);
                fre_all[i] = frequency_matrix (sequence_aa_1, 5, aa[i]->length+extend_len, aa[i]->closure_rows);
                fre2_all[i] = frequency_matrix (sequence_aa_2, 5, 2*(aa[i]->length+extend_len), aa[i]->closure_rows);
                sco_all[i] = get_profile (sequence_aa_1, 5, aa[i]->length+extend_len, aa[i]->closure_rows);
                /*improve the profile by decreasing the middle part of motif*/
                /*sco_all[i] = impovre_profle (sco_all[i],aa[i]->length);*/

                sco2_all[i] = get_profile (sequence_aa_2, 5, 2*(aa[i]->length+extend_len), aa[i]->closure_rows);
        }
        double length_local_1;
	int min_length = 0;
        continuous *I_11, *I_22;
        for (ii=0,IC_1=0;ii<(closure_id-1);ii++,IC_1=0)
        {
                if (IS_duplicate[ii]) continue;
                length_local_1 = aa[ii]->length+extend_len;

                I_11 = get_column_IC (aa, ii);
                for (j=floor((aa[ii]->length+extend_len)/2); j<floor((aa[ii]->length+extend_len)/2)+length_local_1;j++)
                {
                        IC_1 = IC_1 + I_11[j]*I_11[j];
                }
                for (jj=ii+1,IC_2=0;jj<closure_id; jj++,IC_2=0,similarity = -10)
                {
                        if (IS_duplicate[ii]) continue;
			min_length = MIN(aa[jj]->length+extend_len, aa[ii]->length+extend_len);
                	length_local_1 = aa[jj]->length+extend_len;
			I_22 = get_column_IC (aa, jj);
                        for (j=floor((aa[jj]->length/2+extend_len)); j<floor((aa[jj]->length+extend_len)/2)+length_local_1;j++)
                        {
                                IC_2 = IC_2 + I_22[j]*I_22[j];
                        }
                        int startpos[2], start, mos=0;
                        continuous *I_12, *I_21;
                        for(startpos[0]=0,startpos[1]=aa[jj]->length+extend_len; startpos[1]>0||startpos[0]<aa[ii]->length+extend_len; mos=0)
                        {
                                I_12 = get_column_IC_12 (aa, ii, jj, startpos[0], startpos[1]);
                                I_21 = get_column_IC_21 (aa, ii, jj, startpos[0], startpos[1]);
                                IC_21=IC_12=0;
                                for (start = mos; start < mos+min_length; start++)
                                {
                                        IC_12 = IC_12 + I_12[start]*I_22[startpos[1]+start];
                                        IC_21 = IC_21 + I_21[start]*I_11[startpos[0]+start];
                                }
				IC_temp = (IC_12 + IC_21)/(IC_1 + IC_2);
                                if (IC_temp > similarity)
                                {
                                        similarity = IC_temp;
                                }
                                for(p=startpos[0]+1,q=startpos[1]+1, IC_21=IC_12=0;p< aa[ii]->length+1+extend_len && q< aa[jj]->length+1+extend_len;p++,q++,IC_21=IC_12=0)
                                {
                                        mos++;
                                        for (start = mos; start < mos+min_length; start++)
                                        {
                                                IC_12 = IC_12 + I_12[start]*I_22[startpos[1]+start];
                                                IC_21 = IC_21 + I_21[start]*I_11[startpos[0]+start];
                                        }
                                        IC_temp = (IC_12 + IC_21)/(IC_1 + IC_2);
                                        if (IC_temp > similarity)
                                        {
                                                similarity = IC_temp;
                                        }
                                }
                                if(startpos[1]>0) startpos[1]--;
                                else startpos[0]++;
                                free (I_12);
                                free (I_21);
                        }
                        free (I_22);
			if (similarity > threshold)
			{
				IS_duplicate[jj] = TRUE;
			}
                }
                free (I_11);
                if ((ii+1)%10==0) {
					uglyTime("cleanup", ii);
					verboseDot();
				}
        }
	for (i=0;i<closure_id;i++)
        {
                for (j=0;j<5;j++)
                {
                        free (fre_all[i][j]);
                        free (sco_all[i][j]);
                }
                free (fre_all[i]);
                free (sco_all[i]);
        }
        free (fre_all);
        free (sco_all);
        for (i=0;i<closure_id;i++)
        {
                for (j=0;j<5;j++)
                {
                        free (fre2_all[i][j]);
                        free (sco2_all[i][j]);
                }
                free (fre2_all[i]);
                free (sco2_all[i]);
        }
        free (fre2_all);
        free (sco2_all);
        return IS_duplicate;
}
/************************************************************************/
discrete **get_closure_matrix_1 (Closures **aa, int closure_id, continuous threshold)
/*in order to reduce time complexity*/
{
	
	bool *IS_duplicate; /*remove duplicate from the closures which are corresponding to a common operon*/
	AllocArray (IS_duplicate,closure_id);
	clo_matr1 = alloc2d (closure_id,closure_id);
	discrete **matrix;/*store the 0-1 matrix used to find clique*/
        matrix = alloc2d (closure_id,closure_id);
	continuous **matrix_score;
	matrix_score = alloc2dd (closure_id,closure_id);
	int i,j,ii,jj;

	for (i=0;i<closure_id;i++)
		IS_duplicate[i] = FALSE;

	int p_opt=0,q_opt=0;
	for (ii=0;ii<(closure_id);ii++)
        {
                for (jj=0;jj<closure_id; jj++)
		{
                        matrix[ii][jj] = 0;
			matrix_score[ii][jj] = 0;
			clo_matr1[ii][jj] = 0;
		}
                matrix[ii][ii]=1;
        }
	int p,q;
	continuous IC_1=0,IC_2=0,IC_12=0,IC_21=0,IC_temp =0,similarity=-1;
	fre_all = (discrete***)malloc(sizeof(discrete**)*closure_id);
	fre2_all = (discrete***)malloc(sizeof(discrete**)*closure_id);
	sco_all = (continuous***)malloc(sizeof(continuous**)*closure_id);
	sco2_all = (continuous***)malloc(sizeof(continuous**)*closure_id);
	char **sequence_aa_1, **sequence_aa_2;
	for (i=0;i<closure_id;i++)
	{
		fre_all[i] = alloc2d (5, aa[i]->length+extend_len);
		fre2_all[i] = alloc2d (5, 2*aa[i]->length+extend_len);
		sco_all[i] = alloc2dd (5, aa[i]->length+extend_len);
		sco2_all[i] = alloc2dd (5, 2*aa[i]->length+extend_len);
	}
	for (i=0;i<closure_id;i++)
	{
		sequence_aa_1 = get_sequences_from_closures (aa[i]);
		sequence_aa_2 = get_2L_sequeces_from_closures (aa[i]);
		fre_all[i] = frequency_matrix (sequence_aa_1, 5, aa[i]->length+extend_len, aa[i]->closure_rows);
		fre2_all[i] = frequency_matrix (sequence_aa_2, 5, 2*(aa[i]->length+extend_len), aa[i]->closure_rows);
		sco_all[i] = get_profile (sequence_aa_1, 5, aa[i]->length+extend_len, aa[i]->closure_rows);
		/*improve the profile by decreasing the middle part of motif*/
		sco_all[i] = impovre_profle (sco_all[i],aa[i]->length+extend_len);

		sco2_all[i] = get_profile (sequence_aa_2, 5, 2*(aa[i]->length+extend_len), aa[i]->closure_rows);
/*		sco2_all[i] = impovre_profle (sco2_all[i],2*aa[i]->length);*/
	}
	int forepart, endpart;
	double length_local_1;
	
	for (ii=0;ii<(closure_id-1);ii++)
        {
                /*printf ("%d\t%s\n",ii+1,aa[ii]->name);*/
		if (IS_duplicate[ii]) continue;
		length_local_1 = aa[ii]->length+extend_len;
		forepart=floor(length_local_1/3),endpart=ceil(length_local_1/3);
		for (i=1;i<5;i++)
                        for (j=0;j<aa[ii]->length+extend_len;j++)
			{
                                IC_1 = IC_1 + ((fre_all[ii][i][j]+p_markov[i][0])*sco_all[ii][i][j])/(aa[ii]->closure_rows+1);
/*				if (ii==3)	printf ("1111111111111\t%d\t%d\t%3.2f\t%3.2f\t%3.2f\n",j,fre_all[ii][i][j],sco_all[ii][i][j],IC_1,IC_2);*/
			}
		for (jj=ii+1;jj<closure_id; jj++)
                {
			/* printf ("%d\t%d\t%3.2f\t%3.2f\n",ii,jj,IC_1,IC_2);*/
			for (i=1;i<5;i++)
                		for (j=0;j<aa[jj]->length+extend_len;j++)
				{	
					IC_2 = IC_2 + ((fre_all[jj][i][j]+p_markov[i][0])*sco_all[jj][i][j])/(aa[jj]->closure_rows+1);
/*					if (ii==3 && jj==8)      printf ("2222222222222222\t%d\t%d\t%3.2f\t%3.2f\t%3.2f\n",j,fre_all[jj][i][j],sco_all[ii][i][j],IC_2,IC_1);*/
				}
/*			printf ("%d\t%d\t%3.2f\t%3.2f\n",ii,jj,IC_1,IC_2);*/
			for (p=0;p<aa[ii]->length+extend_len+1;p++)
			{
		                for (q=0; q< aa[jj]->length+extend_len+1; q++)
		                {
                		        for (i=1;i<5;i++)
					{
/*						for (j=0;j<length_local_1;j++)
                                                {
                                                        IC_12 = IC_12 + ((fre2_all[ii][i][p+j]+p_markov[i][0])*sco2_all[jj][i][q+j])/(aa[ii]->closure_rows+1);
                                                        IC_21 = IC_21 + ((fre2_all[jj][i][q+j]+p_markov[i][0])*sco2_all[ii][i][p+j])/(aa[jj]->closure_rows+1);
                                                }*/
						/* we improve the profile by descrease the middle part of motif*/
				                for (j=0;j<forepart;j++)
				                {
                        				IC_12 = IC_12 + ((fre2_all[ii][i][p+j]+p_markov[i][0])*sco2_all[jj][i][q+j])/(aa[ii]->closure_rows+1);
							IC_21 = IC_21 + ((fre2_all[jj][i][q+j]+p_markov[i][0])*sco2_all[ii][i][p+j])/(aa[jj]->closure_rows+1);
				                }
						for (j=forepart;j<length_local_1-endpart;j++)
                                                {
                                                        IC_12 = IC_12 + ((fre2_all[ii][i][p+j]+p_markov[i][0])*sco2_all[jj][i][q+j])/(2*aa[ii]->closure_rows+1);
                                                        IC_21 = IC_21 + ((fre2_all[jj][i][q+j]+p_markov[i][0])*sco2_all[ii][i][p+j])/(2*aa[jj]->closure_rows+1);
                                                }
						for (j=length_local_1-endpart; j<length_local_1; j++)
						{
							IC_12 = IC_12 + ((fre2_all[ii][i][p+j]+p_markov[i][0])*sco2_all[jj][i][q+j])/(aa[ii]->closure_rows+1);
                                                        IC_21 = IC_21 + ((fre2_all[jj][i][q+j]+p_markov[i][0])*sco2_all[ii][i][p+j])/(aa[jj]->closure_rows+1);
						}
					}
                		        IC_temp = (IC_12+IC_21)/(IC_1+IC_2);
/*					if (ii==0 && jj==5) printf ("%d\t%d\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n",p,q,IC_1,IC_2,IC_12,IC_21,IC_temp);*/
		                        if (IC_temp > similarity)  
					{
						similarity = IC_temp;
						p_opt = p;
						q_opt = q;
					}
					IC_12 = IC_21 =0;
		                }
			}
/*			printf ("%d\t%s\t%d\t%s\t%3.2f\n",ii,aa[ii]->name,jj,aa[jj]->name,similarity);*/
			if (po->IS_reference)
			{
				/*we clean up the closures which have a similarity score over 0.8 when vertical*/
				if ((similarity>0.5) && (sameString(aa[ii]->name,aa[jj]->name ))) 
				{
					IS_duplicate[jj] = TRUE;
					IC_2 =0;
					continue;
				}
				if (similarity < 0) 
				{
					IC_2 =0;
					continue;
				}
				/*we increase the weight of similarity score across different corresponding operons*/
				if (!sameString(aa[ii]->name,aa[jj]->name)) similarity = 1.2*similarity; 
			}
			/*update the matrix and clo_matr1 which save 0-1 and relative position, respectively*/
/*			printf ("%d\t%s\t%d\t%s\t%3.2f\n",ii,aa[ii]->name,jj,aa[jj]->name,similarity);*/
			matrix_score[ii][jj] = similarity;
			matrix_score[jj][ii] = similarity;
			clo_matr1[ii][jj] = p_opt - q_opt;
	                clo_matr1[jj][ii] = q_opt - p_opt;
/*			printf ("%d\t%d\t%d\n",ii,jj,clo_matr1[ii][jj]);*/

			/*when we combine the closures get from horizonal, we use threshold to filter similarity*/
			if (po->IS_reference_H)
			{
				if (similarity > threshold)
        	                {
                		        matrix[ii][jj]=1;
                                	matrix[jj][ii]=1;
			/*		printf ("%d\t%d\t%3.2f\n",ii,jj,similarity);*/
					clo_matr1[ii][jj] = p_opt - q_opt;
					clo_matr1[jj][ii] = q_opt - p_opt;
                        	}
			}
			/*if there are no -H and -g, we use threshold to filter similarity*/
			if (!po->IS_reference_H && !po->IS_reference)
			{
				if (similarity > threshold)
                                {
                                        matrix[ii][jj]=1;
                                        matrix[jj][ii]=1;
/*                                        printf ("%d\t%d\t%3.2f\n",ii,jj,similarity);*/
                                        clo_matr1[ii][jj] = p_opt - q_opt;
                                        clo_matr1[jj][ii] = q_opt - p_opt;
                                }
			}
			IC_2=0;
			similarity = -1;
                }
		IC_1=0;
		if ((ii+1)%10==0) 
		{
			uglyTime("get closure %d", ii);
			verboseDot();
		}
        }
	/*print out the similarity scores*/
/*	printf ("\no\t");
	for (i=0;i<closure_id;i++)
		printf ("%d\t",i);
	printf ("\n");
	for (i=0;i<closure_id;i++)
	{
		printf ("%d\t",i);
		for (j=0;j<closure_id;j++)
			printf ("%3.2f\t",matrix_score[i][j]);
		printf ("\n");
	}*/
	/* when combine the closures from footprinting data, update the 0-1 matrix without a threshold*/
	if (po->IS_reference) update_matrix (matrix_score,matrix,closure_id);
	/*print out the matrix*/
/*        printf ("\no\t");
        for (i=0;i<closure_id;i++)
                printf ("%d\t",i);
        printf ("\n");
        for (i=0;i<closure_id;i++)
        {
                printf ("%d\t",i);
                for (j=0;j<closure_id;j++)
                        printf ("%d\t",matrix[i][j]);
                printf ("\n");
        }*/
	
	/*print out the matrix*/
 /*       printf ("\no\t");
        for (i=0;i<closure_id;i++)
                printf ("%d\t",i);
        printf ("\n");
        for (i=0;i<closure_id;i++)
        {
                printf ("%d\t",i);
                for (j=0;j<closure_id;j++)
   		        printf ("%d\t",clo_matr1[i][j]);
                printf ("\n");
        }*/

	for (i=0;i<closure_id;i++)
	{
		for (j=0;j<5;j++)
		{
			free (fre_all[i][j]);
			free (sco_all[i][j]);
		}
		free (fre_all[i]);
                free (sco_all[i]);
	}
	free (fre_all);
        free (sco_all);
	for (i=0;i<closure_id;i++)
        {
                for (j=0;j<5;j++)
                {
                        free (fre2_all[i][j]);
                        free (sco2_all[i][j]);
                }
                free (fre2_all[i]);
                free (sco2_all[i]);
        }
        free (fre2_all);
        free (sco2_all);
	free (IS_duplicate);
        return matrix;
}
/************************************************************************/
void update_matrix (continuous **matrix_score,discrete **matrix,int num)
/*update matrix base on matrix_score*/
{
	int i,j,num_1,num_2,jj,max,pos=0,num_a,num_b;
	int top = MIN (floor(num/2)+1, 40);
/*	printf ("top is %d\n",top);*/
	for (i=0;i<num;i++)
	{
		num_a = get_genome_num_from_closure (all[i]);
		num_1 = 0;
		for (j=0;j<num;j++)
		{
			if (matrix[i][j] == 1) 
			{
				num_1++;
/*				printf ("1111111111%d\t%d\t%3.2f\n",i,j,matrix_score[i][j]);*/
			}
		}
		if (num_1 >= top) continue;
/*		printf ("**%d\t%d\n",num_1,i);*/
		for (j=0;j<num;j++)
		{
/*						printf ("%d\t%d\t%d\t%d\n",i,j,genome[i]->seq_num,genome[j]->seq_num);*/
			if (matrix_score[i][j]>=0.85)
			{
				num_b = get_genome_num_from_closure (all[j]);
			/*	printf ("%d\t%d\t%d\t%d\n",i,j,num_a,num_b);*/
				matrix[i][j] = 1;
				matrix[j][i] = 1;
				matrix_score[i][j] =0;
/*				printf ("22222222222%d\t%d\t%3.2f\n",i,j,matrix_score[i][j]);*/
				matrix_score[j][i] =0;
				num_1++;
			}
		}
		if (num_1 >= top) continue;
		for (j=0;j<num;j++)
		{
			if (matrix_score[i][j] <= 0.5)
			{
				matrix_score[i][j] = 0;
				matrix_score[j][i] = 0;
			}
		}
		for (j=num_1;j<top;j++)
		{
			max = 0;
			num_2 = 0;
			for (jj=0;jj<num;jj++)	
			{
				if (matrix_score[i][jj] == 0) 
				{
					num_2++;
					continue;
				}
				if (matrix_score[i][jj] > max)
				{
					max = matrix_score[i][jj];
					pos = jj;
				}
			}
/*			printf ("%d\t%d\t%d\n",pos,num_2,num);*/
			if (num_2 == num) break;
			num_b = get_genome_num_from_closure (all[pos]);
/*					printf ("%d\t%d\t%d\t%d\n",i,pos,num_a,num_b);*/
			if (is_two_closures_regulon (all[i], all[pos], genome[num_a], genome[num_b], i, pos))
			{
			/*printf ("%d\t%d\t%d\t%d\n",i,j,num_a,num_b);*/
				matrix[i][pos] = 1;
				matrix[pos][i] = 1;
				matrix_score[i][pos] =0;
				matrix_score[pos][i] =0;
/*			printf ("3333333333333333%d\t%d\t%d\t%3.2f\n",i,j,pos,matrix_score[i][pos]);*/
			}
			matrix_score[i][pos] =0;
			matrix_score[pos][i] =0;
		}
	}
}
/************************************************************************/
discrete **get_closure_matrix (Closures **aa, int closure_id, continuous threshold)
/*get the 0-1 matrix base on similarity of closures so that we can find clique of closures*/
{
	discrete **matrix;
	matrix = alloc2d (closure_id,closure_id);
	int ii, jj;
	for (ii=0;ii<(closure_id);ii++)
	{
		for (jj=0;jj<closure_id; jj++)
			matrix[ii][jj]=0;
		matrix[ii][ii]=1;
	}
	continuous similarity;
	for (ii=0;ii<(closure_id-1);ii++)
        	for (jj=ii+1;jj<closure_id; jj++)
                {
			similarity =  similarity_closures(aa[ii], aa[jj]);
			if (similarity > threshold)
			{
				matrix[ii][jj]=1;
				matrix[jj][ii]=1;
			}
                }
	return matrix;
}
/************************************************************************/
void closure_clique (const char* fn)
{
        FILE *fw = mustOpen(fn, "w");
        int i, j, cnt,cnt_1=0;
        int rec_num = 0;
        bool *ver_c;
	ver = clo_num;
	arr_c = alloc2d (ver,ver);
	for (i=0;i<ver;i++)
		for (j=0;j<ver;j++)
			arr_c[i][j] = clo_matr[i][j];
        AllocArray(ver_c, ver);
        for (i=0;i<ver;i++)
            ver_c[i]=TRUE;
        for (i=0;i<ver;i++)
        {
                cnt=0;
                for (j=0;j<ver;j++)
                        cnt+=arr_c[i][j];
                if (cnt < po->COL_WIDTH-1)
                {
                        ver_c[i]=FALSE;
			cnt_1++;
                }
        }
        /*if there is not enough information, print the error message and exit */
        if (cnt_1 == ver)
        {
                printf ("\nSorry, the generated closures can not be combined any more\n");
                exit(1);
        }
        /* edge_ptr describe edges */
        AllocArray(edge_list, HEAP_SIZE);
        /* Allocating heap structure */
        struct fibheap *heap;
        heap = fh_makeheap();
        fh_setcmp(heap, edge_cmpr);
        /* Generating seed list and push into heap */
        Edge __cur_min = {0, 0, po->COL_WIDTH-1};
        Edge *_cur_min = &__cur_min;
        Edge **cur_min = & _cur_min;
	for (i = 0; i < ver; i++)
	{
		if (ver_c[i])
                        for (j = i+1; j < ver; j++)
                                if (ver_c[j])
                                {
                                        cnt = str_intersect_r(arr_c[i], arr_c[j]);
					if (cnt < (*cur_min)->score)  continue;
                                        AllocVar(edge_ptr);
                                        edge_ptr -> gene_one = i;
                                        edge_ptr -> gene_two = j;
                                        edge_ptr -> score = cnt;
                                        fh_insert_fixed(heap, edge_ptr, cur_min);
                                }
	}
	rec_num = heap->fh_n;
	if (rec_num == 0)
	{
		printf ("\nSorry, the generated closures can not be combined any more\n");
                exit(1);
	}
        /* sort the seeds */
        ReAllocArray(edge_list, rec_num);
        fh_dump(heap, edge_list);
        /* motif seeds finding (clustering)*/
        int n_regulon = 0;
        progress("\nRegulon finding (refinement and expansion) started");
        n_regulon = regulon(fw, edge_list,rec_num);
        uglyTime("\n%d Regulons are written to %s", n_regulon, fn);
        /* clean up */
        for (i=0; i<rec_num; i++)
                free(edge_list[i]);
        free(edge_list);
	free(ver_c);
}

/************************************************************************/
int regulon (FILE *fw, Edge **el, int n)
{
        AllocArray(genes, clo_num);
        AllocArray(conds, clo_num);
        int i;
        for (i=0;i<clo_num;i++)
        {
                genes[i]=i;
                conds[i]=i;
        }
	int block_id = 0;
        Block **bb;
        int allocated = po->SCH_BLOCK;
        AllocArray(bb, allocated);
        Edge *e;
        Block *b;
        struct dyStack *genes, *scores, *b_genes, *allincluster;
        int  j, k, components;
        AllocArray(profile, clo_num);
        for (j = 0; j < clo_num; j++) 
		AllocArray(profile[j], 2);
        genes = dsNew(clo_num);
        scores = dsNew(clo_num);
        allincluster = dsNew(clo_num);
        bool *candidates;
        AllocArray(candidates, clo_num);

        e = *el; i = 0;
        while (i++ < n)
        {
                e = *el++;
                /* check if any one of the two genes already enumerated in previous blocks */
                bool flag = TRUE;
                if ( isInStack(allincluster,e->gene_one)|| isInStack(allincluster,e->gene_two)||arr_c[e->gene_one][e->gene_two]==0) flag = FALSE;
                if (!flag) continue;
                for (j = 0; j < clo_num; j++)
                        for (k = 0; k < 2; k++) 
				profile[j][k] = 0;
                AllocVar(b);
                b->score = MIN(2, e->score);
                /* initialize the stacks genes and scores */
                int ii;
                dsClear(genes);
                dsClear(scores);
                for(ii = 0; ii < clo_num; ii ++)
                {
                        dsPush(genes,-1);
                        dsPush(scores,-1);
                }
                dsClear(genes);
                dsClear(scores);
                dsPush(genes, e->gene_one);
                dsPush(genes, e->gene_two);
                dsPush(scores, 1);
                dsPush(scores, b->score);
		/* branch-and-cut condition for seed expansion */
                int cand_threshold = floor(po->COL_WIDTH * po->TOLERANCE);
                if (cand_threshold < 2) cand_threshold = 2;
                /* maintain a candidate list to avoid looping through all rows */
                for (j = 0; j < clo_num; j++) 
			candidates[j] = TRUE;
                candidates[e->gene_one] = candidates[e->gene_two] = FALSE;
                components = 2;
                /* expansion step, generate a bicluster without noise */
                block_init(e, b, genes, scores, candidates, cand_threshold, &components, allincluster);
                /* track back to find the best score that which genes makes it */
                for(k = 0; k < components; k++)
                        if ((dsItem(scores,k) == b->score)&&(dsItem(scores,k+1)!= b->score)) break;
                components = k + 1;
                int ki;
                for (ki=0; ki < clo_num; ki++)
                	candidates[ki] = TRUE;
                for (ki=0; ki < components - 1 ; ki++)
                {
                        seed_update(arr_c[dsItem(genes,ki)]);
                        candidates[dsItem(genes,ki)] = FALSE;
                }
                candidates[dsItem(genes,k)] = FALSE;
                genes->top = k ;
                int cnt = 0;
                bool *colcand;
                AllocArray(colcand, clo_num);
                for(ki = 0; ki < clo_num; ki++) 
			colcand[ki] = FALSE;
                /* add columns satisfy the conservative r */
                seed_current_modify(arr_c[dsItem(genes,k)], colcand, &cnt, components);
                /* add some new possible genes */
                int m_cnt;
                for ( ki = 0; ki < rows; ki++)
                {
			/* the element appeared in allincluster can be used again if there is -r*/
                        if (!po->IS_closure)
				if (isInStack(allincluster,ki)) candidates[ki] = FALSE;
                        m_cnt = intersect_row(colcand, arr_c[dsItem(genes,0)], arr_c[ki]);
                        if ( candidates[ki] && (m_cnt >= floor(cnt*po->TOLERANCE)) )
                        {
                                dsPush(genes,ki);
                                components++;
                                candidates[ki] = FALSE;
                        }
                }
                b->block_rows_pre = components;
                free(colcand);
                /* save the current cluster*/
		b_genes = dsNew(b->block_rows_pre);
                for (ki = 0; ki < b->block_rows_pre; ki++)
                        dsPush(b_genes, dsItem(genes,ki));
                /* store gene arrays inside block */
		b->genes = dsNew(components);
                b->conds = dsNew(clo_num);
                scan_block(b_genes, b);
                if (b->block_cols == 0) continue;
                b->block_rows = components;
                b->score = b->block_rows;
                dsClear(b->genes);
                for ( ki=0; ki < components; ki++)
                        dsPush(b->genes,dsItem(genes,ki));
                for(ki = 0; ki < components; ki++)
                        if(!isInStack(allincluster, dsItem(genes,ki))) dsPush(allincluster,dsItem(genes,ki));
		/* count the operon number of predicted regulons if po->IS_reference is TRUE*/
		bool *reference;
		int oper_num = 0;
		if (po->IS_reference_H) b->oper_num = 1;
		AllocArray (reference,components);
		for (ki=0; ki<components; ki++) reference[ki] = FALSE;
		if (po->IS_reference)
		{
			for (ki=0;ki<components; ki++)
			{
				if (!reference[ki]) 
				{
					reference[ki]= TRUE; 
					for (k=ki+1;k<components;k++)
					{
						if (strcmp(all[dsItem(b->genes,ki)]->name, all[dsItem(b->genes,k)]->name)==0)
							reference[k]= TRUE;

					}
					oper_num++;
				}
			}
			b->oper_num = oper_num;
		/*	printf ("%d\t%d\n",b->oper_num,b->score);*/
		}
		/* save current block to global **bb */
                bb[block_id++] = b;
                /* reaching the results number limit */
                if (block_id == po->SCH_BLOCK) break;
                /*verboseDot();*/
        }
        /* free-up the candidate list */
        free(candidates);
        free(allincluster);
        return report_regulon(fw, bb, block_id);
}

/************************************************************************/
int get_num_TF (Closures *aa)
{
	int jj,ii=0,kk;
	if (!po->IS_SWITCH) return 1;
	for (jj=0; jj<aa->closure_rows; jj++)
	{
		int i1=0, i2=0;
		for (kk = 0; kk<anno[dsItem(aa->sequence,jj)]->num;kk++)
                {
                	i1 = dsItem(anno[dsItem(aa->sequence,jj)]->init,kk);
                        i2 = dsItem(anno[dsItem(aa->sequence,jj)]->end,kk);
                        if ((dsItem(aa->position,jj) >= i1)&&(dsItem(aa->position,jj) <= i2))  	ii++;        
		}
	}
	return ii;
}
/************************************************************************/
void print_frequency_anno (FILE *fw, char **aa, int num, int num1)
/*print the STA of TFBS frequency*/
{
	fprintf (fw,"----------------------------------------------------\nSTA_%d: ", num1);
	int i,j;
	bool *isanno;
	AllocArray (isanno,num);
	for (i=0;i<num;i++)
		isanno[i]=FALSE;
	for (i=0;i<num;i++)
	{
		if (isanno[i]) continue;
		isanno[i]= TRUE;
		int k=1;
		for (j=i+1;j<num;j++)
			if (sameString(aa[i],aa[j])&&(!isanno[j]))
			{	
				k++;
				isanno[j] = TRUE;
			}
		/*we do not output the regulons annotation whose TFBSs is less than 5*/
		if (k > 2) fprintf (fw,"%s_%d ",aa[i],k);
	}
	fprintf (fw, "\n");
}
/************************************************************************/
int get_genome_num_from_closure (Closures *aa)
{
	int i, genome_num=0;
	bool flag;
	for (i=0;i<reference_oper_num;i++)
        {
                if (sameString(genome[i]->oper_name,aa->name))
                {
                        flag = TRUE;
                        genome_num = i;
                        break;
                }
        }
	return genome_num;
}
/************************************************************************/
