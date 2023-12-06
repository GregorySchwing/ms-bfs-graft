/*
 Pothen-Fan Algorithm: DFS-based maximum cardinality matching algorithm in bipartite graphs.
 Developed by Ariful Azad
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <list>
#include <vector>
#include <omp.h>
#include <stdlib.h>
#include "graphgenBP.h"
#include "maximalMatching.h"
using namespace std;


//#define FULL_STAT 1 // uncomment it if you want more stats


// ---------- maximum matching routines -------------------
int* Pothen_Fan_Fairnes(graph* G, int* mateI);



int main(int argc, char** argv)
{
	if(argc != 3)
	{
		printf("Usage: ./pf filename numThreads\n");
		return -1;
	}
	
	
	int numThreads = atoi(argv[2]);
    omp_set_num_threads(numThreads);
	char inFile[200];
	strcpy(inFile,argv[1]); 
	
	// Create graph
	graph* g = (graph *) malloc(sizeof(graph));
	process_mtx_compressed(inFile, g);
    
    
    
    int isolated_rows = 0, isolated_cols = 0;
#pragma omp parallel
    {
        int tisor = 0, tisoc = 0; //thread private variables
#pragma omp for
        for(int u=0; u<g->nrows; u++)
        if(g->vtx_pointer[u+1] == g->vtx_pointer[u]) tisor++;
#pragma omp for
        for(int u=g->nrows; u<g->n; u++)
        if(g->vtx_pointer[u+1] == g->vtx_pointer[u]) tisoc++;
        
        __sync_fetch_and_add(&isolated_rows,tisor);
        __sync_fetch_and_add(&isolated_cols,tisoc);
    }
    
    
    printf("\n===================================\n");
    printf("Problem Statistics\n");
    printf("===================================\n");
    printf("vertices : %ld\n", g->n);
    printf("rows = %ld cols = %ld\n", g->nrows, g->n - g->nrows);
    printf("Isolated rows = %ld (%.2lf %%)\n", isolated_rows, (double)100 * isolated_rows/g->nrows);
    printf("Isolated cols = %ld (%.2lf %%)\n", isolated_cols, (double)100* isolated_cols/(g->n - g->nrows));
    printf("Number of edges : %ld\n", g->m);
    printf("Number of threads: %d\n", numThreads);
    printf("===================================\n");


	
	
	int NV = g-> n;
	int nrows = g-> nrows;
	int* unmatchedU = (int*) malloc(NV * sizeof(int));
	int* mateI = (int*) malloc(NV * sizeof(int));
	for(int u=0; u<NV; u++)
	{
		mateI[u] = -1;
	}
	int numUnmatchedU;
	numUnmatchedU = KarpSipserInitS(g, unmatchedU,  mateI);
    //int* mate = Pothen_Fan_Fairnes(g, mateI);
    //free(mate);
	
    
    int threads[]={1,2,4,8,15,30,60,120,240};
    int* mate;
    for(int i=0; i<9; i++)
    {
        omp_set_num_threads(threads[i]);
        mate = Pothen_Fan_Fairnes(g, mateI);
        free (mate);
    }
	
	free_graph(g);
	free(g);
	return 0;
}





// DFS with lookahead that finds a single augmenting path 
// called from pothen-fan
int findAugPathLookahead(int sFirst, int* flagLookahead, int* flagDFS, int* mate, graph* G,int* path,int* edgeIndexLookahead, int* edgeIndex, int* edgeVisited)
{
	
	int *edgeStart = G->vtx_pointer;
	int *endVertex = G->endV;
	int NE = G->m;
	int top = -1;
	path[++top] = sFirst; // push , path is equivalent to stack
	
	while (top >= 0 )// while stack not empty
	{
		int u = path[top];
		int uDegree = edgeStart[u+1] - edgeStart[u];
		// lookahed part
		while(++edgeIndexLookahead[u] < uDegree)
		{
			int v = endVertex[edgeStart[u]+ edgeIndexLookahead[u]];
			(*edgeVisited) ++;  // just for stat
			if(__sync_fetch_and_add(&flagLookahead[v],1) == 0)
			{
				if(mate[v] == -1)
				{
					__sync_fetch_and_add(&flagDFS[v],1);
					path[++top] = v; // push
					return top+1; // top = augmenting path length
					
				}
				
			}
		}
		
		while(++edgeIndex[u] < uDegree)
		{
			
			int v = endVertex[edgeStart[u]+ edgeIndex[u]];
			(*edgeVisited) ++;
			if(__sync_fetch_and_add(&flagDFS[v],1) == 0)
			{
				if(mate[v] != -1) // means other vertex already allocate this in lookahed phase
				{
					
					path[++top] = v; // push v
					path[++top] = mate[v]; // push next u
					break;
					
				}
			}
			
		}
		if(edgeIndex[u]==uDegree)
		{
			top-= 2;// pop
		}
	}
	return top+1;
}





// DFS with lookahead that finds a single augmenting path 
// called from pothen-fan
int findAugPathLookaheadReverse(int sFirst, int* flagLookahead, int* flagDFS, int* mate, graph* G,int* path,int* edgeIndexLookahead, int* edgeIndex, int* edgeVisited)
{
	
	int *edgeStart = G->vtx_pointer;
	int *endVertex = G->endV;
	int NE = G->m;
	int top = -1;
	path[++top] = sFirst; // push , path is equivalent to stack 
	
	while (top >= 0 )// while stack not empty 
	{
		int u = path[top];
		int uDegree = edgeStart[u+1] - edgeStart[u];
		// lookahed part
		while(++edgeIndexLookahead[u] < uDegree)
		{
			int v = endVertex[edgeStart[u]+ edgeIndexLookahead[u]];
			(*edgeVisited) ++;  // just for stat
			if(__sync_fetch_and_add(&flagLookahead[v],1) == 0)
			{
				if(mate[v] == -1)
				{
					__sync_fetch_and_add(&flagDFS[v],1);
					path[++top] = v; // push
					return top+1; // top = augmenting path length
					
				}
				
			}
		}
		
		//while(++edgeIndex[u] < uDegree)
		while(--edgeIndex[u] >= 0)
		{
			
			int v = endVertex[edgeStart[u]+ edgeIndex[u]];
			(*edgeVisited) ++;
			if(__sync_fetch_and_add(&flagDFS[v],1) == 0) 
			{
				if(mate[v] != -1) // means other vertex already allocate this in lookahed phase 
				{
					
					path[++top] = v; // push v
					path[++top] = mate[v]; // push next u
					break;
					
				}
			}
			
		}
		if(edgeIndex[u]==-1)
		{
			top-= 2;// pop
		}
		
		
	}
	return top+1;
}

// ------------- PF with Fairness ---------------------

int* Pothen_Fan_Fairnes(graph* G, int* mateI)
{
	
	double time2,time;
	//time = omp_get_wtime();
	int NE = G->m;
	int NV = G->n;
	int *endVertex = G->endV;
	int *edgeStart = G->vtx_pointer;
	
	
	int* unmatchedU = (int*) malloc(NV * sizeof(int));
	int* tQ = (int*) malloc(NV * sizeof(int));
	int* mate = (int*) malloc(NV * sizeof(int));
	int* flagDFS = (int*) malloc(NV * sizeof(int));
	int* flagLookahead = (int*) malloc(NV * sizeof(int));
	int* edgeIndex = (int*) malloc(NV * sizeof(int));
	int* edgeIndexLookahead = (int*) malloc(NV * sizeof(int));
	time = omp_get_wtime();
    vector<vector< double> > fullStat;
	
    
#define THREAD_BUF_LEN 16384
    int numUnmatchedU = 0;
    
    // identify unmatched and non-isolated vertices from where search will begin
#pragma omp parallel
    {
        int kbuf=0, nbuf[THREAD_BUF_LEN];
#pragma omp for
        for(int u=0; u<G->nrows; u++)
        {
            if(mateI[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
            {
                if (kbuf < THREAD_BUF_LEN)
                {
                    nbuf[kbuf++] = u;
                }
                else
                {
                    int voff = __sync_fetch_and_add (&numUnmatchedU, THREAD_BUF_LEN);
                    for (int vk = 0; vk < THREAD_BUF_LEN; ++vk)
                    unmatchedU[voff + vk] = nbuf[vk];
                    nbuf[0] = u;
                    kbuf = 1;
                }
            }
        }
        if(kbuf>0)
        {
            int voff = __sync_fetch_and_add (&numUnmatchedU, kbuf);
            for (int vk = 0; vk < kbuf; ++vk)
            unmatchedU[voff + vk] = nbuf[vk];
        }
    }
    
    
#pragma omp parallel for default(shared)
	for(int i=0; i<NV; i++)
	{
		mate[i] = mateI[i];
		flagLookahead[i] = 0;
		edgeIndexLookahead[i] = -1;      
	}
	
	
	int nthreads;
#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}
	int** augPaths = (int**) malloc(sizeof(int*) * nthreads);
	for(int i=0; i<nthreads; i++)
	{
		augPaths[i] = (int*) malloc(NV * sizeof(int));
	}
	
	
	int iterations = 0;
	int numEdgeVisited = 0;
	
    printf("\n************* Starting Pothen-Fan Algorithm  *************\n");
    printf("Initial number of non-isolated row vertices = %ld\n\n", numUnmatchedU);

    printf(" ====================Phase by phase statistics===============================\n");
    printf(" Phase   Initial-unmatched  Matched-in-this-phase    MaxAugPathLen  Time (sec)\n");
    printf(" ============================================================================\n");
    
    
	while(1)
	{
		iterations++;
		time2 = omp_get_wtime();
		int tQ_len = 0;
		if(iterations % 2 == 1) // odd iterations
		{
#pragma omp parallel for schedule(static)
			for(int i=0; i< NV; i++)
			{
				flagDFS[i] = 0;
				edgeIndex[i] = -1;
			}
		}
		else
		{
#pragma omp parallel for schedule(static)
			for(int i=0; i< NV; i++)
			{
				flagDFS[i] = 0;
				edgeIndex[i] = edgeStart[i+1] - edgeStart[i];
			}
		}
		int maxAugPathLen = -1;
		int phaseEdgeVisited = 0;
#pragma omp parallel for schedule(dynamic) default(shared)
		for(int i=0; i < numUnmatchedU; i++)
		{
			
#ifdef FULL_STAT
			double timeTree = omp_get_wtime();
#endif
			int tid = omp_get_thread_num();
			int* augPath = augPaths[tid];
			int edgeVisited = 0;
			int uFirst = unmatchedU[i];
			int augPathLen;
			if(iterations % 2 == 1) // odd iterations
				augPathLen = findAugPathLookahead(uFirst,flagLookahead, flagDFS,mate,G,augPath,edgeIndexLookahead,edgeIndex, &edgeVisited) ;
			else
				augPathLen = findAugPathLookaheadReverse(uFirst,flagLookahead, flagDFS,mate,G,augPath,edgeIndexLookahead,edgeIndex, &edgeVisited) ;
			if (augPathLen > 0)
			{
				// augment in serial ... can be done in parallel also ...
				int u = unmatchedU[i];
				for(int k=0; k< augPathLen; k+=2)
				{
					mate[augPath[k]] = augPath[k+1];
					mate[augPath[k+1]] = augPath[k];
				}
				
			}
			else
			{
				tQ[__sync_fetch_and_add(&tQ_len,1)] = uFirst;
			}
			
            if(augPathLen > maxAugPathLen) maxAugPathLen = augPathLen;
#ifdef FULL_STAT
			__sync_fetch_and_add(&phaseEdgeVisited,edgeVisited);
			
#endif
		}
		numEdgeVisited+= phaseEdgeVisited;
		
		double dfsTime=omp_get_wtime() - time2;
        
        printf("%4ld    %12ld %20ld %18ld %14.2lf\n",iterations, numUnmatchedU, numUnmatchedU - tQ_len, maxAugPathLen, dfsTime);
        
#ifdef FULL_STAT
        double curStat[] = {numUnmatchedU, numUnmatchedU - tQ_len, maxAugPathLen, dfsTime};
        std::vector<double> temp (curStat, curStat + sizeof(curStat) / sizeof(double) );
        fullStat.push_back(temp);
#endif

        if( (tQ_len ==0) || (numUnmatchedU == tQ_len))
        {
            numUnmatchedU  = 0; // every non-isolated rows are matched, just for correct stats
            break;
        }
		int* tt = tQ;
		tQ =  unmatchedU;
		unmatchedU = tt;
		
		numUnmatchedU = tQ_len;
		
	}
	
	double totalTime =omp_get_wtime() - time;
	//printf("DFS with lookahed (PF) Iterations=%d numEdgeVisited=%d Quality=%lf%% total time = %f \n\n",iterations, 100.0*(G->nrows-numUnmatchedU)/G->nrows, numEdgeVisited, totalTime);
   
    // numUnmatchedU contains only non-isolated unmatched vertices
    // compute actual matching cardinality
    int matched_rows = 0;
#pragma omp parallel
    {
        int tmatched = 0; //thread private variables
#pragma omp for
        for(int u=0; u<G->nrows; u++)
        if(mate[u]!=-1) tmatched++;
        __sync_fetch_and_add(&matched_rows,tmatched);
    }
    
    int isolated_rows = G->nrows - matched_rows - numUnmatchedU;
    
    printf("============================================================================\n\n");
    printf("========= Overall Statistics ===============\n");
    printf("Number of  Iterations           = %ld \n", iterations);
    //printf("Avg. Length of Augmenting Paths = %.2lf \n", (double)total_aug_path_len/total_aug_path_count);
    printf("Total time                      = %.2lf sec\n", totalTime);
    printf("Maximum matching cardinality    = %ld (%.2lf%%)\n", matched_rows*2, (double)100.0*(matched_rows*2)/G->n);
    printf("Matched Rows cardinality        = %ld (%.2lf%%)\n", matched_rows, (double)100.0*(matched_rows)/G->nrows);
    printf("Isolated Rows                   = %ld (%.2lf%%)\n", isolated_rows, (double)100.0*(isolated_rows)/G->nrows);
    
    printf("===========================================\n");

    

#ifdef FULL_STAT
    FILE* fp1 = fopen("fullStat.txt","w");
    for (int i=0; i<fullStat.size(); i++)
    {
        fprintf(fp1,"%8ld %8ld %15ld %15ld    %6.4lf\n", (int)fullStat[i][0], (int)fullStat[i][1], (int)fullStat[i][2], (int)fullStat[i][3], fullStat[i][4]);
    }
    fclose(fp1);
#endif


	free(flagDFS);
	free(edgeIndex);
	free(unmatchedU);
	free(flagLookahead);
	free(edgeIndexLookahead);
	free(tQ);
	for(int i=0; i<nthreads; i++)
	{
		free(augPaths[i]);
	}
	free(augPaths);
    
    return (mate);
}






