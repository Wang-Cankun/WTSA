/************************************************************************/
/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Feb. 15, 2010
 * Blocks finding procedure base on the input microarray data
 */

#include "matrix_process.h"

/**************************************************************************/
int post_processing_blocks_1 ( FILE *fw1, Block** bb, int num)
{
        sort_block_list(bb, num);
        int i, j,k;
        int n = MIN(num, po->RPT_BLOCK);
        bool flag;
        Block **output;
        AllocArray(output, n);
        Block **bb_ptr = output;
        Block *b_ptr;
        int cur_rows, cur_cols;
        int inter_rows, inter_cols;
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
                        print_bc_1( fw1, b_ptr, j++);
                        *bb_ptr++ = b_ptr;
                        verboseDot();
                }
        }
        return j;
}
/**************************************************************************/
void print_bc_1 (FILE *fw1, Block* b, int num)
{
	int i, j;
        int block_rows, block_cols;

        /* block height (genes) */
        block_rows = b->block_rows;
        /* block_width (conditions) */
        block_cols = b->block_cols;
        fprintf(fw1, "BC%03d\tS=%d\n", num, block_rows * block_cols);

        fprintf(fw1, " Genes [%d]: ", block_rows);
        for (i=0; i<dsSize(b->genes); i++)
                fprintf(fw1, "%s ", genes_n[dsItem(b->genes, i)]);
        fprintf(fw1, "\n");

        fprintf(fw1, " Conds [%d]: ", block_cols);
        for (i=0; i<dsSize(b->conds); i++)
                fprintf(fw1, "%s ", conds_n[dsItem(b->conds, i)]);
        fprintf(fw1, "\n");
        /* the complete block data output */
        for (i=0; i<dsSize(b->genes); i++)
        {
                fprintf(fw1,"%10s:",genes_n[dsItem(b->genes, i)]);
                for (j=0; j<dsSize(b->conds); j++)
                {
                        fprintf(fw1, "\t%d", symbols[arr_c[dsItem(b->genes, i)][dsItem(b->conds, j)]]);
                }
                fputc('\n', fw1);
                if (i == b->block_rows_pre -1)
                        fputc('\n',fw1);
        }
        fputc('\n', fw1);

}
/**************************************************************************/
int cluster_1 (FILE *fw1, Edge **el, int n)
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
        return post_processing_blocks_1(fw1, bb, block_id);
}
/**************************************************************************/
