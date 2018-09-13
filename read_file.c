/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 16, 2010
 * read all kinds of input files	
 */

#include "read_file.h"
#define MAXC 100000
#define MAX_PROMOTER_NUM 200
#define MAX_NAME_NUM 2000

static int bb[USHRT_MAX];
static char delims[] = "\t\r\n>";
static char delims1[] = ":";
static char *atom = NULL;
/***********************************************************************/
/* read the annotation file, if any, with the format of NC_000913_motif */
void read_annotation(FILE* fp)
{
        char *buffer;
        int i=0;
        buffer=(char*)malloc(sizeof(char)*MAXC);
        /* get the sequences number*/
        while(fgets(buffer,MAXC,fp)!=NULL)  i++;
        Annotation *aa;
        /* make point fp to the begining of the file by rewind function */
	rewind(fp);
        i=0;
        while(fgets(buffer,MAXC,fp)!=NULL)
        {
                int j=0,k=0;
                AllocVar(aa);
                aa->init = dsNew(20);
                aa->end = dsNew(20);
                aa->TF = alloc2c (20,10);
                atom = strtok(buffer, delims);
                atom = strtok(NULL, delims);
                atom = strtok(NULL, delims);
                atom = strtok(NULL, delims);
                while (atom != NULL)
                {
                        j++;
                        if (j%3 == 1) strcpy(aa->TF[k], atom);
                        else if (j%3 == 2) dsPush(aa->init,atoi(atom));
                        else if (j%3 == 0) {dsPush(aa->end, atoi(atom)); k++;}
                        atom = strtok(NULL, delims);
                }
		aa->num = k;
		/* save current annotations to global structure **anno */
                anno[i] = aa;
                i++;
        }
        free (buffer);
	fseek(fp, 0, 0);
}

/***********************************************************************/
/*read the .closures file output by BoBro*/
void read_closures (FILE* fc)
{
        char *buffer;
	double significance=0;
        int i=0,j;
        /*get the number of closures*/
        clo_num=0;
        buffer=(char*)malloc(sizeof(char)*MAXC);
        while(fgets(buffer,MAXC,fc)!=NULL)
                if (buffer[0] == '#')  clo_num++;
        /* we need minus the first five head lines with # at the begining in .closures file */
	clo_num = clo_num-5;
	/* pop up the real number of candidate closures in input .closures file */
        printf ("There are %d candidate closures\n",clo_num);

        AllocArray(all, clo_num);
        /* set the indicator to the begining of the file and skip the first five header rows*/
        rewind(fc);
        for (i=0;i<5;i++)
		fgets(buffer,MAXC,fc);
        discrete *clo_row;
        AllocArray(clo_row, clo_num);
        /*count the number of motif (TFBS) in each closure which are saved in *clo_row*/
        i=0;
        int i_1=0;
        while(fgets(buffer,MAXC,fc)!=NULL)
        {
                if (buffer[0] == '#')
                {
                        j=0;
                        while(fgets(buffer,MAXC,fc)!=NULL)
                        {
                                if (buffer[0] == '-') break;
                                else  j++;
                        }
                        clo_row[i] = j;
                        i++;
                }
		if (buffer[0] == 'I')  i_1++;
        }
        /*read the information of each closure*/
        rewind(fc);
        for (i=0;i<5;i++)
                fgets(buffer,MAXC,fc);
        i=0;
	char *oper_name;
	AllocArray (oper_name,MAX_NAME_NUM);
        Closures *cctemp;
        while(fgets(buffer,MAXC,fc)!=NULL)
        {
		if (buffer[0] == 'I' )
		{
			AllocArray (oper_name,MAX_NAME_NUM);
			atom = strtok(buffer, delims);
			atom = strtok(NULL, delims);
			atom = strtok(NULL, delims);
			strcpy(oper_name, atom);	
		}
		/*read in the pvalue of each closure*/
		if (buffer[0] == ' ' && buffer[1] == 'P')
		{
			atom = strtok(buffer, delims);
			atom = strtok(NULL, delims);
			significance = atof(atom);
		}
		if (buffer[0] == '#')
                {
                        AllocVar(cctemp);
                        int temp = clo_row[i];
                        cctemp->sequence = dsNew(temp);
                        cctemp->position = dsNew(temp);
                        j=0;
			while(fgets(buffer,MAXC,fc)!=NULL)
                        {
                                if (buffer[0] == '-') break;
                                else
                                {
                                        atom = strtok(buffer, delims);
                                        cctemp->length = po->MOTIFLENGTH;
                                        dsPush(cctemp->sequence,atoi(atom));
                                        atom = strtok(NULL, delims);
                                        dsPush(cctemp->position,atoi(atom));
                                        j++;
                                }
                        }
			AllocArray (cctemp->name,MAX_NAME_NUM);
			strcpy(cctemp->name,oper_name);
			cctemp->significance = significance;
/*			printf ("%3.2e\n",cctemp->significance);*/
                        cctemp->closure_rows = j;
			/* save current closure's information to global **all */
                        all[i] = cctemp;
                        i++;
                }
        }
	/* compare similarity between the readed closures and generate a matrix clo_matr 
	 * so that we can find clique base on the 0-1 matrix by function 'regulon' */
        clo_matr = get_closure_matrix_1 (all, clo_num, po->thre);
	fseek(fc, 0, 0);
}
/***********************************************************************/
discrete *change_AGCT_to_num (char *seq, int length)
{
	int i=0;
	discrete *number;
	AllocArray (number,length);
	for (i=0;i<length;i++)
		number[i]=0;
	for (i=0;i<length;i++)
	{
		if (seq[i]=='G' || seq[i]=='g') number[i] = 1;
		if (seq[i]=='C' || seq[i]=='c') number[i] = 2;
		if (seq[i]=='T' || seq[i]=='t') number[i] = 3;
	}
	return number;
}
/***********************************************************************/
int scan_genome (continuous **scoreM, continuous AveScore, int motif_length, FILE* fp)
/*scan each seed in background genome*/
{
		
        sum_genome = 0;
        int background_num =0;
        int i,j,seq_length=0;
        bool IS_end_sequences = FALSE;
        continuous score_scan, score_scan_p;
        continuous **scoreM_panlindromic;
        scoreM_panlindromic = get_palindromic_profile (scoreM, motif_length) ;
        char *buffer, *buffer_combine_end,  *end_sequence;
        buffer=(char*)malloc(sizeof(char)*MAX_SEQUENCE_LENGTH);
        buffer_combine_end =(char*)malloc(sizeof(char)*MAX_SEQUENCE_LENGTH);
        end_sequence =(char*)malloc(sizeof(char)*motif_length);
	discrete *buffer_combine_end_number;
	buffer_combine_end_number = change_AGCT_to_num(buffer_combine_end,seq_length);

        rewind(fp);
        while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp)!=NULL)
        {
                if (buffer[0]=='>') continue;
                if (IS_end_sequences)
                {
                        j=0;
                        for (i=0;i<strlen(end_sequence);i++)
                                buffer_combine_end[j++] = end_sequence[i];
                        for (i=0;i<strlen(buffer)-1;i++)
                        {
                                buffer_combine_end[j++] = buffer[i];
                                sum_genome++;
                        }
                        seq_length = j;
                        /*printf ("all\t%Zu\t%Zu\t%d\n",strlen(buffer_combine_end),strlen(buffer),seq_length);*/
                }
                else
                {
                        j=0;
                        for (i=0;i<strlen(buffer)-1;i++)
                        {
                                buffer_combine_end[j++] = buffer[i];
                                sum_genome++;
                        }
			seq_length = j;
                }
		buffer_combine_end_number = change_AGCT_to_num(buffer_combine_end,seq_length);
                for (i=0;i<seq_length-motif_length+1;i++)
                {
                        score_scan = 0;
                        score_scan_p = 0;
                        for (j=0;j<motif_length;j++)
                        {
				score_scan += scoreM[buffer_combine_end_number[i+j]+1][j];
				score_scan_p += scoreM_panlindromic[buffer_combine_end_number[i+j]+1][j];
                        }
                        if (score_scan >= po->conserve_background*AveScore || score_scan_p >= po->conserve_background*AveScore)
                                background_num++;
                }
                j=0;
                for (i=seq_length-motif_length+1; i<seq_length; i++)
                {
                        end_sequence[j++] = buffer_combine_end[i];
                        IS_end_sequences = TRUE;
                }
	}
        free (buffer);
        free (buffer_combine_end);
        free (end_sequence);
	free (buffer_combine_end_number);
	uglyTime("scan_genome", s_rows);
        return background_num;
}

/***********************************************************************/
void read_sequences(FILE* fp1 )
/*read fasta sequence file*/
{
        char *buffer;
        int i=0,j=0,k=0,t=0;
	/*for check palindromic*/
		char *Palin;
		Palin=(char*)malloc(sizeof(char)*po->PalinLength);
	
	buffer=(char*)malloc(sizeof(char)*MAX_SEQUENCE_LENGTH);
        /* get the sequences number*/
        s_rows=0;
        while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp1)!=NULL)
                if (buffer[0]=='>') s_rows++;
	/* if there are less than two sequences in input file, err pop up */
	if (s_rows < 3)
	{
	       printf ("\nSorry, there is not enough sequences in your fasta file\n");
               exit(1);	   
	}    
	/* save the locus_id*/
	if (po->ID)
	{
		locus_id = alloc2c(s_rows,20);
		rewind(fp1);
		i=0;
		while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp1)!=NULL)
	        {
			if (buffer[0]=='>')
			{
				atom = strtok(buffer, delims);
				atom = strtok(NULL, delims);
        	                atom = strtok(NULL, delims);
				strcpy(locus_id[i++], atom);
			}
		}
	}
	
	/*store the sequence information*/
	SequenceInfo = alloc2c(s_rows,MAX_NAME_NUM);
	rewind(fp1);
        int SequenceId=0;
        while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp1)!=NULL)
        {
                if (buffer[0]=='>')
                {
                        atom = strtok(buffer, delims);
                        if (strlen(atom)>2000)
                        errAbort("\nSorry, the title of squence is too long (>2000 chars)\n");
                        strcpy(SequenceInfo[SequenceId++], atom);
                }
        }

	/*printf ("%d\n",s_rows);*/
	/*read the weighted sequences*/
	if (po->SequenceWeight)
        {
                AllocArray (SequenceWeight, s_rows);
                rewind(fp1);
                i=0;
                while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp1)!=NULL)
                {
                        if (buffer[0]=='>')
                        {
                                atom = strtok(buffer, delims1);
                                atom = strtok(NULL, delims1);
                                SequenceWeight[i++] = atof(atom);
                        }
                }
        }

        s_col=(int*)malloc(sizeof(int)*s_rows);
        promoter_length=(int*)malloc(sizeof(int)*s_rows);
        /* read in the sequences*/
        rewind(fp1);
        sequences= alloc2c(s_rows,MAX_SEQUENCE_LENGTH);
	height_matrix= alloc2d(s_rows,MAX_SEQUENCE_LENGTH);
	
	/* change A T C G to 0 1 2 3*/
        seq_matrix = alloc2d(s_rows,MAX_SEQUENCE_LENGTH);

        /*int temp = strlen(seq_matrix[0]);*/
	for (i=0;i<s_rows; i++)
        {
                for (j=0;j<MAX_SEQUENCE_LENGTH; j++)
                {
                        seq_matrix[i][j] = 0;
                        height_matrix[i][j] = 0;
                }
        }
		
			


        while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp1)!=NULL)
        {
		/* case insensitive */
                char s[1] = ",";
                char *token =NULL;
		if(buffer[0]=='A'||buffer[0]=='T'||buffer[0]=='G'||buffer[0]=='C'||buffer[0]=='a'||buffer[0]=='t'||buffer[0]=='g'||buffer[0]=='c'||buffer[0]=='N'||buffer[0]=='n')
                {
                        t=1;
			/* fix the problem that there may exist \n between two '>' */
                        for (j=k;j<strlen(buffer)-1+k;j++)
			{
				sequences[i][j] = buffer[j-k];
				if (buffer[j-k]=='T' || buffer[j-k]=='t') seq_matrix[i][j] = 3;
				if (buffer[j-k]=='C' || buffer[j-k]=='c') seq_matrix[i][j] = 2;
				if (buffer[j-k]=='G' || buffer[j-k]=='g') seq_matrix[i][j] = 1;
			}
                        k=k+strlen(buffer)-1;
                }
                else if(buffer[0]>= '0' && buffer[0] <= '9' )
                {
                        /*printf("%d\n",strlen(buffer));
                        token = strtok(buffer, s);*/
                        token = strtok(buffer, s);
                        j=0;
                        while (token&&j<k+strlen(buffer)-1) {
                        if(atoi(token)==0) {
                        height_matrix[i][j] = 1;
                        }
                        else {
                                height_matrix[i][j] = atoi(token);
                                
                        }
                        
                        /*printf("%s\t\n",token);
                        printf("%d\t%d\t%d\t\n",i,j,height_matrix[i][j]);*/
                        token = strtok(NULL, ",");
                        j++;
                        }
                        
                        /*for (j=0;j<strlen(buffer);j++)
			{       
                                /*
                                if(i==0)
                                {   
                                token = strtok(NULL,s);
                                }
                                else
                                {
                                height[i][j] = atoi(token);
                                
                                printf("%d\n",atoi(token));
                                
				
			}*/
                }
                else
                {
                        if (t==0) i=0;
                        else i++;
                        k=0;
                }
        }
        s_cols=strlen(sequences[0]);
        /*consider the situation of comparing sequences with different length*/
        for (i=0;i<s_rows;i++)
        {
                s_col[i]=strlen(sequences[i]);
                promoter_length[i]=strlen(sequences[i]);
		if(po->checkPalin)
		{
			/*printf ("4: %s\n", sequences[i]);*/
			if(s_col[i] < po->PalinLength)
			{
				/*printf ("%s: NA\n", SequenceInfo[i]);*/
				continue;
			}
			for (j=0; j<s_col[i]-po->PalinLength+1; j++)
			{
				t=0;
				for (k=j; k<j+po->PalinLength; k++)
				{
					Palin[t] = sequences[i][k];
					t++;
				}
				if (palindromic_pro_1 (Palin, po->PalinLength))
				/*printf ("5: %d-%c-%d-%d: %.2f\n%s\n",j, sequences[i][j],po->PalinLength,t,PalinPro,Palin);*/
					printf ("%s\t%d\t%d\t%s\n",SequenceInfo[i],j+1,po->PalinLength,Palin);
			}
		}
		if (s_col[i] < po->MOTIFLENGTH && !po->checkPalin)
		{
			printf ("\nThe %dth sequence are too short to find motif\n", i+1);
			/*printf ("\nSorry, the %dth sequence are too short to find motif\n", i+1);
	                exit(1);*/
		}
		/*printf ("%d\n",s_col[i]);*/
                if (strlen(sequences[i])>s_cols) 
		{
			s_cols=strlen(sequences[i]);
		}
        }
	if(po->checkPalin)
	{
		exit(1);
	}

	/*generate the fre_matrix*/
	fre_matrix = (discrete***)malloc(sizeof(discrete**)*s_rows);
	for (i=0;i<s_rows;i++)
	{
		fre_matrix[i] = alloc2d (4, s_cols);
		for (j=0;j<4;j++)
			for (k=0;k<s_cols;k++)
				fre_matrix[i][j][k] = 0;
                for (k=0;k<(s_cols-8);k++)
                {
			fre_matrix[i][seq_matrix[i][k]][k] = 1;
                }
	}

        /*the markov processing, generate a markov matrix so that we can simulate sequences*/
        p_markov = markov(sequences,s_rows);
        td_markov = d_markov(sequences,s_rows);
        free (buffer);
	fseek(fp1, 0, 0);
}


/***********************************************************************/
/*read a combined file containing reference genomes for all the confirmed operons*/
void read_reference_genome (FILE* fp)
{
	char *buffer;
	buffer=(char*)malloc(sizeof(char)*MAX_SEQUENCE_LENGTH);
	int oper_num = 0;
	int  promo_num=0, t, k,i,j;
	/* count the number of operons */
	while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp)!=NULL)
                if ((buffer[0]=='>')&&(buffer[1]=='>')) oper_num++;
	/* we minus the last >>end from the counter */	
	reference_oper_num = oper_num-1;
	AllocArray (genome,(oper_num-1));
	oper_num_all = oper_num-1;
	Reference_genome *genome_temp;
	AllocVar (genome_temp);
	oper_num=0;

	char **sequences_temp;
	int *clo_length;
	char name[MAX_NAME_NUM];
	sequences_temp = alloc2c(MAX_PROMOTER_NUM , MAX_SEQUENCE_LENGTH);
	AllocArray (clo_length, MAX_PROMOTER_NUM);
	rewind(fp);
	t=0,k=0;
	while(fgets(buffer,MAX_SEQUENCE_LENGTH,fp)!=NULL)
	{
		if ((buffer[0]=='>')&&(buffer[1]=='>')&& (t==0))
		{
			atom = strtok(buffer, delims);
  			strcpy(name, atom);
			AllocVar (genome_temp);
			AllocArray (genome_temp->oper_name, MAX_NAME_NUM);
			promo_num = 0;
			continue;
		}
		if ((buffer[0]=='>')&&(buffer[1]=='>')&&(t==1))
		{
			genome_temp->sequences_r = alloc2c ((promo_num+1), MAX_SEQUENCE_LENGTH);
			AllocArray (genome_temp->clo_num,(promo_num+1));
			for (i=0; i<promo_num+1; i++)
				genome_temp->clo_num[i] = clo_length [i];
			for (i=0; i<promo_num+1; i++)
				for (j=0; j<clo_length[promo_num]; j++)
					genome_temp->sequences_r[i][j] = sequences_temp[i][j];
			strcpy (genome_temp->oper_name, name);
		       	genome_temp->seq_num = promo_num+1;	
			genome_temp->markov_matrix = markov(genome_temp->sequences_r,genome_temp->seq_num);
			/* save current reference genome to global **genome */
			genome[oper_num++] = genome_temp;
			AllocVar (genome_temp);
			AllocArray (genome_temp->oper_name, MAX_NAME_NUM);
			atom = strtok(buffer, delims);
			strcpy(name, atom);
			promo_num = 0;
			t=0;
			continue;
		}
		/* case insensitive and adopt \n in each sequence */
		if(buffer[0]=='A'||buffer[0]=='T'||buffer[0]=='G'||buffer[0]=='C'||buffer[0]=='a'||buffer[0]=='t'||buffer[0]=='g'||buffer[0]=='c')
                {
                        t=1;
/*			printf ("%Zu\n",strlen(buffer));*/
                        for (j=k;j<strlen(buffer)-1+k;j++)
                                sequences_temp[promo_num][j] = buffer[j-k];
                        k=k+strlen(buffer)-1;
			clo_length [promo_num] = strlen(sequences_temp[promo_num]);
                }
                else
                {
                        if (t==0) promo_num = 0;
                        else promo_num++;
                        k=0;
                }
	}
	fseek(fp, 0, 0);
}
/***********************************************************************/
/*the following functions are used to read microarray data matrix*/
/* Pre-read the datafile, retrieve gene labels and condition labels
 *  as well as determine the matrix size */
void get_matrix_size (FILE* fp)
{
        size_t n = 0;
        char *line;
        if (getline(&line, &n, fp)>=0)
        {
                atom = strtok(line, delims);
                atom = strtok(NULL, delims);
                while (atom != NULL)
                {
                        atom = strtok(NULL, delims);
                        cols++;
                }
        }
        while (getline(&line, &n, fp)>=0)
        {
                atom = strtok(line, delims);
                rows++;
        }
        fseek(fp, 0, 0);
}
/* Read in the labels on x and y, in microarray terms, genes(rows) and conditions(cols)*/
void read_labels (FILE* fp)
{
        int row = 0;
        int col;
        size_t n = 0;
        char *line;
        while (getline(&line, &n, fp)>=0)
        {
                atom = strtok(line, delims);
                if (row >= 1)
                        strcpy(genes_n[row-1], atom);
                atom = strtok(NULL, delims);
                col = 0;
                while (atom != NULL)
                {
                        if (row == 0)
                                strcpy(conds_n[col], atom);
                        atom = strtok(NULL, delims);
                        if (++col == cols) break;
                }
                if (++row == rows+1) break;
        }
        fseek(fp, 0, 0);
}

static int charset_add(discrete *ar, discrete s)
{
        int ps = s + SHRT_MAX;
        if (bb[ps]<0)
        {
                bb[ps] = sigma;
                ar[sigma++] = s;
        }
        return bb[ps];
}
/* initialize data for discretization */
void init_dis_m()
{
        int row, col;
        AllocArray(symbols, USHRT_MAX);
        memset(bb, -1, USHRT_MAX*sizeof(*bb));
        charset_add(symbols, 0);
        arr_c = alloc2d (rows,cols);
        for (row = 0; row < rows; row++)
                for (col = 0; col < cols; col++)
                        arr_c[row][col] = 0;
}
/*read discreted microarray data matrix*/
void read_discrete (FILE* fp)
{
        int row, col, i;
        init_dis_m();
        size_t n = 0;
        char *line;
        row = 1;
	/* Skip first line with condition labels */
        getline(&line, &n, fp);
        while (getline(&line, &n, fp)>=0)
        {
                atom = strtok(line, delims);
		/*skip the first column*/
                atom = strtok(NULL, delims);
                col = 0;
                while (atom != NULL)
                {
                        arr_c[row-1][col] = charset_add(symbols, atoi(atom));
                        atom = strtok(NULL, delims);
                        if (++col == cols) break;
                }
                if (++row == rows+1) break;
        }
        printf("Discretized data contains %d classes with charset [ ", sigma);
        for (i=0;i<sigma;i++)
                printf("%d ", symbols[i]); printf("]\n");
        fseek(fp, 0, 0);
}

/***********************************************************************/

