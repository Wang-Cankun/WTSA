/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 16, 2010
 * make the graph G base on the 0-1 matrix got from compare_sequence    
 */

#include "make_graph.h"
/*we can reduce the HEAP_SIZE when the data contain so many genes taht memory is not enough*/
static const int HEAP_SIZE = 30000000;

/**************************************************************************/

/* String intersection function without string copying, only numbers */
/*caculate the weight of the edge in the first graph*/
int str_intersect_r (const discrete *s1, const discrete *s2)
{
	int common_cnt = 0;
	/* s1 and s2 of equal length, so we check s1 only */
        int i;
        for (i=0;i<rows;i++)
	{
		if (*s1==*s2 && (*s1!=0)) common_cnt++;
		s1++; s2++;
	}
	return common_cnt;
}

/**************************************************************************/

/* Fibonacci heap related subroutines */
int edge_cmpr(void *a, void *b)
{
	int score_a, score_b;
	score_a = ((Edge *)a)->score;
	score_b = ((Edge *)b)->score;
	if (score_a < score_b) 
		return -1;
	if (score_a == score_b) 
		return 0;
	return 1;
}
/* Maintain a fixed size heap */
void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min)
{
	if (a->fh_n < HEAP_SIZE) 
	{
		fh_insert(a, (void *)i);
	}
	else
	{
		if (edge_cmpr(cur_min, i) < 0)
		{
			/* Remove least value and renew */
			fh_extractmin(a);
			fh_insert(a, (void *)i);
			/* Keep a memory of the current min */
			*cur_min = (Edge *)fh_min(a);
		}
	}
}
void fh_dump(struct fibheap *a, Edge **res)
{
	int i;
	int n = a->fh_n;
	for (i=n-1; i>=0; i--)
		res[i] = (Edge *) fh_extractmin(a);
}

/**************************************************************************/

/*void make_graph (Closures **cc, const char* fn, const char* fn1)*/
void make_graph (const char* fn1)
{
	FILE *fw1 = mustOpen(fn1, "w");
	int i, j, cnt,cnt_1=0;
	int rec_num = 0;
	bool *ver_c;
	AllocArray(ver_c, ver);
	for (i=0;i<rows;i++)
	    ver_c[i]=TRUE;
	for (i=0;i<rows;i++)
	{
		cnt=0;
		for (j=0;j<cols;j++)
			cnt+=arr_c[i][j];
		if (cnt < po->COL_WIDTH)  
		{
			ver_c[i]=FALSE;
			cnt_1++;
		}
	}
	/*uglyTime("graph1", s_rows);*/
        /*if there is not enough information, print the error message and exit */
	if (cnt_1 == rows) 
        {
		printf ("\nSorry, there is no significant Motifs (legnth=%d) in your fasta file\n", po->MOTIFLENGTH);
		if (all_id > 0)
		{	
			print_params(fw1);
			j = report_closures(fw1,all,all_id,anno);
		}
		exit(1);
        }	
	/*uglyTime("graph2", s_rows);*/
	/* edge_ptr describe edges */
	AllocArray(edge_list, HEAP_SIZE);
        /* Allocating heap structure */
        struct fibheap *heap;
        heap = fh_makeheap();
        fh_setcmp(heap, edge_cmpr);
        /* Generating seed list and push into heap */
        Edge __cur_min = {0, 0, po->COL_WIDTH};
        Edge *_cur_min = &__cur_min;
        Edge **cur_min = & _cur_min;
	/* iterate over all vertices to retrieve all edges */
	for (i = 0; i < rows; i++)
	{
		if (ver_c[i])	
			for (j = i+1; j < cols; j++)
				if (ver_c[j])
				{
					cnt = str_intersect_r(arr_c[i], arr_c[j]);
					if (cnt < (*cur_min)->score) continue;
					AllocVar(edge_ptr);
					edge_ptr -> gene_one = i;
					edge_ptr -> gene_two = j;
					edge_ptr -> score = cnt;
					/*printf ("\n%d\t%d\t%d\n", edge_ptr -> gene_one,edge_ptr -> gene_two,edge_ptr -> score);*/
					fh_insert_fixed(heap, edge_ptr, cur_min);
				}	
	}
	/*uglyTime("graph3", s_rows);*/
        rec_num = heap->fh_n;
	if (rec_num == 0)
	{
        printf ("\nSorry, there is no significant Motifs (length = %d) in your fasta file\n", po->MOTIFLENGTH);
		if (all_id > 0)
		{
			print_params(fw1);
			uglyTime("report2%d", extend_len);
			j = report_closures(fw1,all,all_id,anno);
		}	
                exit(1);
        }

        /* sort the seeds */
        ReAllocArray(edge_list, rec_num);
        fh_dump(heap, edge_list);
        /* motif seeds finding (clustering)*/
        int n_blocks = 0;
	progress("\nMotif finding started");
	n_blocks = cluster(fw1,edge_list,rec_num);
        /*fflush(fw);*/
	

	/* clean up */
	if ((po->no_enhance || po->FastVersion) && (po->MOTIFLENGTH == po->Up)&& is_extended==1)
	{
		printf("clean memory\t\n");
		for (i=0; i<rec_num; i++)
			free(edge_list[i]);
		free(edge_list);
		for (i=0; i<rows; i++)
		{
			free (arr_c[i]);
			free (arr_c1[i]);
		}
		free(arr_c);
		free(arr_c1);
		
	}
	if (po->no_enhance && po->MOTIFLENGTH == po->Up && is_extended==1)	uglyTime("\nMotif length extended to %d\n %d Motifs are written to %s",extend_len+po->MOTIFLENGTH, n_blocks, fn1);
        
}

/***************************************************************************/
