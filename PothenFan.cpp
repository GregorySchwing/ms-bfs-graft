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
long* Pothen_Fan_Fairnes(graph* G, long* mateI);



int main(int argc, char** argv)
{
	if(argc != 3)
	{
		printf("Usage: ./pf filename numThreads\n");
		return -1;
	}
	
	
	long numThreads = atoi(argv[2]);
    omp_set_num_threads(numThreads);
	char inFile[200];
	strcpy(inFile,argv[1]); 
	
	// Create graph
	graph* g = (graph *) malloc(sizeof(graph));
	process_mtx_compressed(inFile, g);
    
    
    
    long isolated_rows = 0, isolated_cols = 0;
#pragma omp parallel
    {
        long tisor = 0, tisoc = 0; //thread private variables
#pragma omp for
        for(long u=0; u<g->nrows; u++)
        if(g->vtx_pointer[u+1] == g->vtx_pointer[u]) tisor++;
#pragma omp for
        for(long u=g->nrows; u<g->n; u++)
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


	
	
	long NV = g-> n;
	long nrows = g-> nrows;
	long* unmatchedU = (long*) malloc(NV * sizeof(long));
	long* mateI = (long*) malloc(NV * sizeof(long));
	for(long u=0; u<NV; u++)
	{
		mateI[u] = -1;
	}
	long numUnmatchedU;
	numUnmatchedU = KarpSipserInitS(g, unmatchedU,  mateI);
    //long* mate = Pothen_Fan_Fairnes(g, mateI);
    //free(mate);
	
    
    int threads[]={1,2,4,8,15,30,60,120,240};
    long* mate;
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
long findAugPathLookahead(long sFirst, long* flagLookahead, long* flagDFS, long* mate, graph* G,long* path,long* edgeIndexLookahead, long* edgeIndex, long* edgeVisited)
{
	
	long *edgeStart = G->vtx_pointer;
	long *endVertex = G->endV;
	long NE = G->m;
	long top = -1;
	path[++top] = sFirst; // push , path is equivalent to stack
	
	while (top >= 0 )// while stack not empty
	{
		long u = path[top];
		long uDegree = edgeStart[u+1] - edgeStart[u];
		// lookahed part
		while(++edgeIndexLookahead[u] < uDegree)
		{
			long v = endVertex[edgeStart[u]+ edgeIndexLookahead[u]];
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
			
			long v = endVertex[edgeStart[u]+ edgeIndex[u]];
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
long findAugPathLookaheadReverse(long sFirst, long* flagLookahead, long* flagDFS, long* mate, graph* G,long* path,long* edgeIndexLookahead, long* edgeIndex, long* edgeVisited)
{
	
	long *edgeStart = G->vtx_pointer;
	long *endVertex = G->endV;
	long NE = G->m;
	long top = -1;
	path[++top] = sFirst; // push , path is equivalent to stack 
	
	while (top >= 0 )// while stack not empty 
	{
		long u = path[top];
		long uDegree = edgeStart[u+1] - edgeStart[u];
		// lookahed part
		while(++edgeIndexLookahead[u] < uDegree)
		{
			long v = endVertex[edgeStart[u]+ edgeIndexLookahead[u]];
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
			
			long v = endVertex[edgeStart[u]+ edgeIndex[u]];
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

long* Pothen_Fan_Fairnes(graph* G, long* mateI)
{
	
	double time2,time;
	//time = omp_get_wtime();
	long NE = G->m;
	long NV = G->n;
	long *endVertex = G->endV;
	long *edgeStart = G->vtx_pointer;
	
	
	long* unmatchedU = (long*) malloc(NV * sizeof(long));
	long* tQ = (long*) malloc(NV * sizeof(long));
	long* mate = (long*) malloc(NV * sizeof(long));
	long* flagDFS = (long*) malloc(NV * sizeof(long));
	long* flagLookahead = (long*) malloc(NV * sizeof(long));
	long* edgeIndex = (long*) malloc(NV * sizeof(long));
	long* edgeIndexLookahead = (long*) malloc(NV * sizeof(long));
	time = omp_get_wtime();
    vector<vector< double> > fullStat;
	
    
#define THREAD_BUF_LEN 16384
    long numUnmatchedU = 0;
    
    // identify unmatched and non-isolated vertices from where search will begin
#pragma omp parallel
    {
        long kbuf=0, nbuf[THREAD_BUF_LEN];
#pragma omp for
        for(long u=0; u<G->nrows; u++)
        {
            if(mateI[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
            {
                if (kbuf < THREAD_BUF_LEN)
                {
                    nbuf[kbuf++] = u;
                }
                else
                {
                    long voff = __sync_fetch_and_add (&numUnmatchedU, THREAD_BUF_LEN);
                    for (long vk = 0; vk < THREAD_BUF_LEN; ++vk)
                    unmatchedU[voff + vk] = nbuf[vk];
                    nbuf[0] = u;
                    kbuf = 1;
                }
            }
        }
        if(kbuf>0)
        {
            long voff = __sync_fetch_and_add (&numUnmatchedU, kbuf);
            for (long vk = 0; vk < kbuf; ++vk)
            unmatchedU[voff + vk] = nbuf[vk];
        }
    }
    
    
#pragma omp parallel for default(shared)
	for(long i=0; i<NV; i++)
	{
		mate[i] = mateI[i];
		flagLookahead[i] = 0;
		edgeIndexLookahead[i] = -1;      
	}
	
	
	long nthreads;
#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}
	long** augPaths = (long**) malloc(sizeof(long*) * nthreads);
	for(long i=0; i<nthreads; i++)
	{
		augPaths[i] = (long*) malloc(NV * sizeof(long));
	}
	
	
	long iterations = 0;
	long numEdgeVisited = 0;
	
    printf("\n************* Starting Pothen-Fan Algorithm  *************\n");
    printf("Initial number of non-isolated row vertices = %ld\n\n", numUnmatchedU);

    printf(" ====================Phase by phase statistics===============================\n");
    printf(" Phase   Initial-unmatched  Matched-in-this-phase    MaxAugPathLen  Time (sec)\n");
    printf(" ============================================================================\n");
    
    
	while(1)
	{
		iterations++;
		time2 = omp_get_wtime();
		long tQ_len = 0;
		if(iterations % 2 == 1) // odd iterations
		{
#pragma omp parallel for schedule(static)
			for(long i=0; i< NV; i++)
			{
				flagDFS[i] = 0;
				edgeIndex[i] = -1;
			}
		}
		else
		{
#pragma omp parallel for schedule(static)
			for(long i=0; i< NV; i++)
			{
				flagDFS[i] = 0;
				edgeIndex[i] = edgeStart[i+1] - edgeStart[i];
			}
		}
		long maxAugPathLen = -1;
		long phaseEdgeVisited = 0;
#pragma omp parallel for schedule(dynamic) default(shared)
		for(long i=0; i < numUnmatchedU; i++)
		{
			
#ifdef FULL_STAT
			double timeTree = omp_get_wtime();
#endif
			long tid = omp_get_thread_num();
			long* augPath = augPaths[tid];
			long edgeVisited = 0;
			long uFirst = unmatchedU[i];
			long augPathLen;
			if(iterations % 2 == 1) // odd iterations
				augPathLen = findAugPathLookahead(uFirst,flagLookahead, flagDFS,mate,G,augPath,edgeIndexLookahead,edgeIndex, &edgeVisited) ;
			else
				augPathLen = findAugPathLookaheadReverse(uFirst,flagLookahead, flagDFS,mate,G,augPath,edgeIndexLookahead,edgeIndex, &edgeVisited) ;
			if (augPathLen > 0)
			{
				// augment in serial ... can be done in parallel also ...
				long u = unmatchedU[i];
				for(long k=0; k< augPathLen; k+=2)
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
		long* tt = tQ;
		tQ =  unmatchedU;
		unmatchedU = tt;
		
		numUnmatchedU = tQ_len;
		
	}
	
	double totalTime =omp_get_wtime() - time;
	//printf("DFS with lookahed (PF) Iterations=%d numEdgeVisited=%d Quality=%lf%% total time = %f \n\n",iterations, 100.0*(G->nrows-numUnmatchedU)/G->nrows, numEdgeVisited, totalTime);
   
    // numUnmatchedU contains only non-isolated unmatched vertices
    // compute actual matching cardinality
    long matched_rows = 0;
#pragma omp parallel
    {
        long tmatched = 0; //thread private variables
#pragma omp for
        for(long u=0; u<G->nrows; u++)
        if(mate[u]!=-1) tmatched++;
        __sync_fetch_and_add(&matched_rows,tmatched);
    }
    
    long isolated_rows = G->nrows - matched_rows - numUnmatchedU;
    
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
    for (long i=0; i<fullStat.size(); i++)
    {
        fprintf(fp1,"%8ld %8ld %15ld %15ld    %6.4lf\n", (long)fullStat[i][0], (long)fullStat[i][1], (long)fullStat[i][2], (long)fullStat[i][3], fullStat[i][4]);
    }
    fclose(fp1);
#endif


	free(flagDFS);
	free(edgeIndex);
	free(unmatchedU);
	free(flagLookahead);
	free(edgeIndexLookahead);
	free(tQ);
	for(long i=0; i<nthreads; i++)
	{
		free(augPaths[i]);
	}
	free(augPaths);
    
    return (mate);
}






