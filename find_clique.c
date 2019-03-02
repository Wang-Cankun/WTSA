/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 15, 2010
 * Clique finding procedure, greedy heuristic by picking an edge with highest
 * score and then dynamically adding vertices into the block and see if 
 * the block score can be improved. 
 * Base on the graph constructed in make_graph.c
 */

#include "find_clique.h"
/************************************************************************/
/* Initialize seed */
static void update_colcand(bool *colcand,struct dyStack *cocluster,discrete *g1, discrete *g2)
{
	int i;
	for (i=0; i< cols; i++)
	{
		if (colcand[i] && (g1[i] != g2[i]))
			colcand[i] = FALSE;
	}
}
/************************************************************************/
int intersect_row(const bool *colcand, discrete *g1, discrete *g2)
{
	int i;
	int cnt = 0;
	for (i=0; i< rows; i++)
		if (colcand[i] && (g1[i] == g2[i]) && (g1[i]!=0)) cnt++;
	return cnt;
}
/************************************************************************/
void seed_update (const discrete *s)
{
        int i;
        for (i = 0; i < cols; i++)
                profile[i][s[i]]++;
}
/************************************************************************/
void seed_current_modify (const discrete *s, bool *colcand, int* cnt, int components)
/* calculate the coverage of any row to the current consensus  cnt = # of valid consensus columns*/
{
	int i, k, flag, n;
	int threshold = components;
	discrete ss;
	*cnt = 0;
	for (i=0; i<cols; i++)
	{   
		flag = 0; ss = s[i];
		for (k=1; k<2; k++)
		{       
			n = profile[i][k];
			if (k == ss) n++;
			if (n >= threshold)
			{
				flag = k; 
				break;
			}
		}
		if (flag)	
		{
			(*cnt)++;
			colcand[i] = TRUE;
		}
	}
}
void block_init(Edge *e, Block *b, 
                     struct dyStack *genes, struct dyStack *scores,
                     bool *candidates, const int cand_threshold,
                     int *components, struct dyStack *allincluster)
{
	int i,score;
	int cnt = 0;
	int max_cnt, max_i;
	struct dyStack *cocluster;
	cocluster=dsNew(cols);
	int *arr_rows, *arr_rows_b;
	AllocArray(arr_rows, rows);
	AllocArray(arr_rows_b, rows);	
	bool *colcand;
	AllocArray(colcand, cols);
	for (i=0; i< cols; i++) 
		colcand[i] = FALSE;
	discrete *g1, *g2;
	g1 = arr_c[dsItem(genes,0)];
	g2 = arr_c[dsItem(genes,1)];
	for (i=0; i< cols; i++)
		if ((g1[i] == g2[i])&&(g1[i]!=0)) 
		{
			colcand[i] = TRUE;
		  	dsPush(cocluster,i);
		}
	for (i = 0; i < rows; i++)
	{
		arr_rows[i] = intersect_row(colcand, arr_c[dsItem(genes,0)], arr_c[i]);
		arr_rows_b[i] = arr_rows[i];
		if (!isInStack(cocluster,i)) candidates[i]= FALSE;
                /*the elements appeared in allin cluster can not be used again when there is no -r*/
		if (!po->IS_closure)
			if (isInStack(allincluster,i)) candidates[i] = FALSE;
	}
	while (*components < rows)
	{
		max_cnt = -1;
		max_i = -1;
		(*components)++;
		for (i=0; i< rows; i++)
		{
			if (!candidates[i]) continue;
			cnt = intersect_row(colcand,arr_c[dsItem(genes,0)],arr_c[i]);
			if (cnt > max_cnt)
			{
				max_cnt = cnt;
				max_i = i;
			}
		}
		if ( max_i < 0) break;
		else
		{	
			score = MIN(*components, max_cnt);
			if (score < b->score) break;
	 		else
			{
				b->score = score;
		 		dsPush(genes, max_i);
		 		dsPush(scores,score);
			 	update_colcand(colcand,cocluster,arr_c[dsItem(genes,0)], arr_c[max_i]);
			 	dsClear(cocluster);
			 	for (i=0; i< cols; i++)
		      			if (colcand[i])	dsPush(cocluster,i);
        	         	candidates[max_i] = FALSE;
                	 	int ii;
		 		for (ii=0;ii<rows;ii++)
			      		if (candidates[ii] && (!isInStack(cocluster,ii))) candidates[ii]= FALSE;
			}
		}
	}
	free(colcand);
	free(cocluster);
}
/******************************************************************/
/* scan through all columns and identify the set within threshold,
 * "fuzziness" of the block is controlled by TOLERANCE (-c)
 */
void scan_block (struct dyStack *gene_set, Block *b_ptr)
{
        int i, j;
        bool flag;
        int block_rows, cur_rows;
        block_rows = cur_rows = dsSize(gene_set);
        int k;
        for (j = 0; j < cols; j++)
                for (k=0; k<2; k++) profile[j][k] = 0;
        for (j = 0; j< cur_rows ; j++)
                seed_update(arr_c[dsItem(gene_set,j)]);
        int btolerance = ceil(po->TOLERANCE*block_rows);
        for (j = 0; j < cols; j++)
        {
                flag = FALSE;
                /* See if this column satisfies tolerance*/
                for (i = 1; i <2; i++)
                        /* Pay attention to the ignored char '.' */
                        if (profile[j][i] >= btolerance)
                        {
                                flag = TRUE; 
				break;
                        }
                if (flag) dsPush(b_ptr->conds, j);
        }
        b_ptr->block_cols = dsSize(b_ptr->conds);
}

/************************************************************************/
/* Core algorithm */
int cluster (FILE *fw1, Edge **el, int n)
{
	
	int block_id = 0;
	Block **bb;
	int allocated = po->SCH_BLOCK;
	AllocArray(bb, allocated);
	Edge *e;
	Block *b;
	struct dyStack *genes, *scores, *b_genes, *allincluster;
	int i, j, k, components;
        AllocArray(profile, cols);
	for (j = 0; j < cols; j++) 
		AllocArray(profile[j], 2/*sigma*/);
	genes = dsNew(rows);
	scores = dsNew(rows);
	allincluster = dsNew(rows);
	bool *candidates;
	AllocArray(candidates, rows);

	e = *el; i = 0;
	while (i++ < n)
	{	
		e = *el++;
		/* check if both genes already enumerated in previous blocks */
		bool flag = TRUE;
		if ( isInStack(allincluster,e->gene_one)|| isInStack(allincluster,e->gene_two)||arr_c[e->gene_one][e->gene_two]==0) flag = FALSE;  
		if (!flag) continue;
		for (j = 0; j < cols; j++)
			for (k = 0; k < 2; k++) 
				profile[j][k] = 0;
		AllocVar(b);
                b->score = MIN(2, e->score);
		/* initialize the stacks genes and scores */		
		int ii;		
		dsClear(genes);
		dsClear(scores);		
		for(ii = 0; ii < rows; ii ++)
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
		for (j = 0; j < rows; j++) 
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
		for (ki=0; ki < rows; ki++)
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
		AllocArray(colcand, cols);
		for(ki = 0; ki < cols; ki++) 
			colcand[ki] = FALSE;             
                /* add columns satisfy the conservative r */ 
		seed_current_modify(arr_c[dsItem(genes,k)], colcand, &cnt, components);
                /* add some new possible genes */
		int m_cnt;
		for ( ki = 0; ki < rows; ki++)
		{
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
		b->conds = dsNew(cols);
		scan_block(b_genes, b);
		if (b->block_cols == 0) continue;
		b->block_rows = components;
                b->score = b->block_rows;
		if (b->score < cand_threshold)
			continue;
		/*printf ("%d\t%d\n",block_id,b->score);*/
		dsClear(b->genes);
		for ( ki=0; ki < components; ki++)
			dsPush(b->genes,dsItem(genes,ki));
		for(ki = 0; ki < components; ki++)
			if(!isInStack(allincluster, dsItem(genes,ki))) dsPush(allincluster,dsItem(genes,ki));	
		bb[block_id++] = b;
                /* reaching the results number limit */
		if (block_id == po->SCH_BLOCK) break;
		/*verboseDot();*/
	}
        /* free-up the candidate list */
	if (po->MOTIFLENGTH == po->Up)
	{
        free(candidates);
		free(allincluster);
	}
	return post_processing_blocks(fw1, bb, block_id);
}


/************************************************************************/
static int post_processing_blocks( FILE *fw1, Block** bb, int num)
{
	/*printf ("%d\n",num);*/
	int closure_id = 0;
	/*uglyTime("post_processing_blocks start %d", closure_id);*/
	sort_block_list(bb, num);
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
			if (inter_rows*inter_cols > po->FILTER*cur_rows*cur_cols)
			{
				flag = FALSE; 
				break;
			}
			k++;
		}
       		i++;
		if (flag)
		{
			print_bc( fw1,cc, closure_id++, b_ptr, j++);
			*bb_ptr++ = b_ptr;
			verboseDot();
		}
	}
	
	if (po->MOTIFLENGTH == po->Up && po->no_enhance)
		print_params(fw1);
	if (po->Low == po->Up)
	{
                for (k=0; k<closure_id;k++)
                {
                        all[all_id] = cc[k];
                        all_id++;
                }
		clo_num = all_id;
		clo_matr = alloc2d (clo_num,clo_num);
		sort_closures_list(all, all_id);
		clo_matr = get_closure_matrix_1 (all, clo_num,po->thre);
		if (po->IS_reference_H)  printf ("\nClosures/Regulons refinement and expansion\n");

		if (po->no_enhance) {
			j = report_closures(fw1,all,all_id,anno);
		}

	}
	else if (po->Low < po->Up)
	{
		int kk = MIN(closure_id,po->number);
		for (k=0; k<kk;k++)
		{
			all[all_id] = cc[k];
			all_id++;
		}
		/*uglyTime("\nWe save top %d Motifs among the found %d Motifs", kk, closure_id);*/
		if (po->MOTIFLENGTH == po->Up && po->no_enhance)
			j = report_closures(fw1,all,all_id,anno);
	}

	return j;
}

/************************************************************************/
int block_cmpr(const void *a, const void *b)
/* compare function for qsort, descending by score */
{
	return ((*(Block **)b)->score - (*(Block **)a)->score);
}

void sort_block_list(Block **el, int n)
{
	qsort(el, n, sizeof *el, block_cmpr);
}
/************************************************************************/
