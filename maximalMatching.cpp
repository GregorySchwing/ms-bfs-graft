/*
 Ariful Azad: Lawrence Berkeley National Laboratory
 */

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "maximalMatching.h"


// helper function used in Serial Karp-Sipser initialization
void findMateS(int u, graph* G, int* flag,int* mate, int* degree)
{
	if(flag[u] != 0) return;
    flag[u] = 1;
	int *endVertex = G->endV;
	int *edgeStart = G->vtx_pointer;
	
	int neighbor_first = edgeStart[u];
	int neighbor_last = edgeStart[u+1];
	for(int j=neighbor_first; j<neighbor_last; j++)
	{
		int v = endVertex[j];
		if(flag[v] == 0) // if I can lock then this v node is unmatched
		{
            flag[v] = 1;
			mate[u] = v;
			mate[v] = u;
			// update degree
			int neighborFirstU = edgeStart[v];
			int neighborLastU = edgeStart[v+1];
			for(int k=neighborFirstU; k< neighborLastU; k++)
			{
				int nextU = endVertex[k];
                degree[nextU]--;
				if( degree[nextU] == 1)
				{
					
					findMateS(nextU,G,flag,mate,degree);
                    
				}
				
			}
			break;
		}
	}
}

// Serial Karp-Sipser maximal matching
int KarpSipserInitS(graph* G, int* unmatchedU,  int* mate)
{
	int nrows = G->nrows;
	int * degree = (int*) malloc(sizeof(int) * nrows);
	int* degree1Vtx = (int*) malloc(sizeof(int) * nrows);
	
	int *endVertex = G->endV;
	int *edgeStart = G->vtx_pointer;
	int nrowsV = G->n;
	int numUnmatchedU = 0;
	int* flag = (int*) malloc(sizeof(int) * nrowsV);
	
	double timeStart = omp_get_wtime();
	
	for(int i=0; i< nrowsV; i++)
	{
		flag[i] = 0;
		mate[i] = -1;
	}
	
	int degree1Tail = 0;
	int degree1Count = 0;
	
	
	
	for(int u=0; u<nrows; u++)
	{
		degree[u] = edgeStart[u+1] - edgeStart[u];
		if(degree[u] == 1)
		{
			degree1Vtx[degree1Count++] = u;
		}
	}
	
	
	
	for(int u=0; u<degree1Count; u++)
	{
            findMateS(degree1Vtx[u],G,flag,mate,degree);
	}
	
	
	for(int u=0; u<nrows; u++)
	{
		if(flag[u] == 0 && degree[u]>0)
			findMateS(u,G,flag,mate,degree);
	}
	
	double timeInit = omp_get_wtime()-timeStart;
	for(int u=0; u<nrows; u++)
	{
		
		if(mate[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
		{
			unmatchedU[numUnmatchedU++] = u;
		}
	}
    
    int matched_rows = 0;
    for(int u=0; u<nrows; u++)
    {
        if(mate[u]!=-1) matched_rows++;
    }

	
    printf("===========================================\n");
    printf("Serial Karp-Sipser Initialization\n");
    printf("===========================================\n");
    printf("Matched Rows        = %ld (%.2lf%%)\n", matched_rows, (double)100.0*(matched_rows)/nrows);
    printf("Computation time    = %lf\n", timeInit);
    printf("===========================================\n");
	//printf("%lf %lf \n", 100.0 * (nrows - numUnmatchedU)/nrows, timeInit);
	free(degree1Vtx);
	free(degree);
	free(flag);
	return numUnmatchedU;
}



// helper function used in Karp-Sipser initialization 
void findMate(int u, graph* G, int* flag,int* mate, int* degree)
{
	if(__sync_fetch_and_add(&flag[u],1) != 0) return;
	int *endVertex = G->endV;
	int *edgeStart = G->vtx_pointer;
	
	int neighbor_first = edgeStart[u];
	int neighbor_last = edgeStart[u+1];   
	for(int j=neighbor_first; j<neighbor_last; j++) 
	{
		int v = endVertex[j];
		if(__sync_fetch_and_add(&flag[v],1) == 0) // if I can lock then this v node is unmatched
		{
			mate[u] = v;
			mate[v] = u;
			// update degree
			int neighborFirstU = edgeStart[v];
			int neighborLastU = edgeStart[v+1];
			for(int k=neighborFirstU; k< neighborLastU; k++)
			{
				int nextU = endVertex[k];
				if( __sync_fetch_and_add(&degree[nextU],-1) == 2)
				{
					
					findMate(nextU,G,flag,mate,degree);
				}
				
			}
			break;
		} 
	}
}


// Multithreaded  Karp-Sipser maximal matching
int KarpSipserInit(graph* G, int* unmatchedU,  int* mate)
{
	int nrows = G->nrows;
	int * degree = (int*) malloc(sizeof(int) * nrows);
	int* degree1Vtx = (int*) malloc(sizeof(int) * nrows);
	
	int *endVertex = G->endV;
	int *edgeStart = G->vtx_pointer;
	int nrowsV = G->n;
	int numUnmatchedU = 0;
	int* flag = (int*) malloc(sizeof(int) * nrowsV);
	
	double timeStart = omp_get_wtime();
	
#pragma omp parallel for default(shared) schedule(static)
	for(int i=0; i< nrowsV; i++)
	{
		flag[i] = 0;
		mate[i] = -1;
	}
	
	int degree1Tail = 0;
	int degree1Count = 0;  
	
	
	
	// populate degree and degree1Vtx
#pragma omp parallel for default(shared) schedule(static)//schedule(dynamic)
	for(int u=0; u<nrows; u++)
	{      
		degree[u] = edgeStart[u+1] - edgeStart[u];
		if(degree[u] == 1)
		{
			degree1Vtx[__sync_fetch_and_add(&degree1Count,1)] = u;
			//flag[u] = 1; // means already taken 
		}
	}
	
	
	
#pragma omp parallel for default(shared) //schedule(dynamic,100)
	for(int u=0; u<degree1Count; u++)
	{
		//findMate1(degree1Vtx[u],G,flag,mate,degree,degree1Vtx,&degree1Head);
		findMate(degree1Vtx[u],G,flag,mate,degree);		  
	}
	
	
	// process other vertices 
#pragma omp parallel for default(shared) schedule(dynamic,100)//schedule(dynamic)
	for(int u=0; u<nrows; u++)
	{
		if(flag[u] == 0 && degree[u]>0)
			findMate(u,G,flag,mate,degree);	  
	}
	
	double timeInit = omp_get_wtime()-timeStart;
#pragma omp parallel for default(shared) 
	for(int u=0; u<nrows; u++)
	{
		
		if(mate[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
		{
			unmatchedU[__sync_fetch_and_add(&numUnmatchedU, 1)] = u;
		}
	}
    
    int matched_rows = 0;
#pragma omp parallel
    {
        int tmatched = 0; //thread private variables
#pragma omp for
        for(int u=0; u<nrows; u++)
        {
            if(mate[u]!=-1) tmatched++;
        }
        __sync_fetch_and_add(&matched_rows,tmatched);
    }

	
    printf("===========================================\n");
    printf("Karp-Sipser Initialization\n");
    printf("===========================================\n");
    printf("Matched Rows        = %ld (%.2lf%%)\n", matched_rows, (double)100.0*(matched_rows)/nrows);
    printf("Computation time    = %lf\n", timeInit);
    //printf("%lf %lf", 100.0 * (nrows - numUnmatchedU)/nrows, timeInit);
    printf("===========================================\n");
	
	free(degree1Vtx);
	free(degree);
	free(flag);
	return numUnmatchedU;
}







