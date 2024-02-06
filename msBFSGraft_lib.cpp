/*=========================================================================
 Multithreaded algorithm for computing maximum Cardinality matching
 in a bipartite graph.
 Author: Ariful Azad (azad@lbl.gov) and Aydin Buluc (abuluc@lbl.gov)
 Please cite: "A Parallel Tree Grafting Algorithm for Maximum Cardinality
 Matching in Bipartite Graphs", A. Azad, A. Buluc, A. Pothen, IPDPS 2015.
 
 *** Copyright Notice ***
 
 "Parallel Maximum Cardinality Matchings via Tree-Grafting" Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
 
 If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
 
 NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
=========================================================================
*/


#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include <algorithm>
#include "graphgenBP.h"
#include "maximalMatching.h"
using namespace std;
#include <sys/time.h>
double getTimeOfDay() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

int* MS_BFS_Graft(graph* G, int* mateI);

extern "C"
int main_lib_msbfsgraft(int argc, char *argv[], int **rows, int **cols, int **matching, int*nr_ptr, int*nc_ptr, int*nn_ptr, int just_read_file, double *parse_graph_time, double *create_csr_time, double *init_time, int * initial_match_count)
{
	if(argc != 4)
	{
		printf("Usage: ./msBFSGraft fileName numThreads parallelKS\n");
		return -1;
	}
    
	int numThreads = atoi(argv[2]);
    int parallelKS = atoi(argv[3]);

    omp_set_num_threads(numThreads);
#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int nprocs   = omp_get_num_procs();
        int tid = omp_get_thread_num();
    }
	
    
    char inFile[100];
    strcpy(inFile,argv[1]);
    FILE* fp;
    fp = fopen(inFile, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error! Could not open input file. Exiting ...\n");
        exit(-1);
    }
    fclose(fp);
    
    graph* g = (graph *) malloc(sizeof(graph));
    //graph* g1 = (graph *) malloc(sizeof(graph));
    double start_parse_graph = getTimeOfDay();
    process_mtx_compressed(inFile, g, rows, cols, matching, nr_ptr, nc_ptr, nn_ptr, parse_graph_time, create_csr_time);
    double end_parse_graph = getTimeOfDay();
    *parse_graph_time = end_parse_graph-start_parse_graph;
    //fast_mtx_read_build(inFile,g);  // ABAB: to replace process_mtx_compressed
    //graph* g =  swap_side(tg);
    //isEqual(g,g1);
    if (just_read_file){
        free(g);
        return 0;
    }
    double start_init = getTimeOfDay();
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

    if (parallelKS){
	    numUnmatchedU = KarpSipserInit(g, unmatchedU,  mateI); //parallel version

    } else {
        numUnmatchedU = KarpSipserInitS(g, unmatchedU,  mateI); //serial version
    }
    
    //int* mate = MS_BFS_Graft(g, mateI); // result is stored in mate array
    *matching = MS_BFS_Graft(g, mateI); // result is stored in mate array
    double end_init = getTimeOfDay();
    *init_time = end_init-start_init;
    // for scaling study
    /*
    int threads[]={1,2,4,8,15,30,60,120,240};
    int* mate;
    for(int i=0; i<9; i++)
    {
        omp_set_num_threads(threads[i]);
        mate = MS_BFS_Graft(g, mateI);
        free (mate);
    }
   */

    #include <vector>
    #include <algorithm>

    // Assuming you have this boolean array to track matched vertices
    std::vector<bool> matched_bool(g->n, false);

    // Declare match_count vector
    std::vector<int> match_count(g->n, 0);

    // Variable to store the total number of matched vertices
    int total_matched_vertices = 0;


    for (int u = 0; u < g->n; u++) {
        // Check for matching based on the conditions u % (g->n/2) and mate[u] % (g->n/2)
        if (mateI[u] > -1 && !matched_bool[u % (g->n/2)] && !matched_bool[mateI[u] % (g->n/2)]) {
            // Update matched information
            matched_bool[u % (g->n/2)] = true;
            matched_bool[mateI[u] % (g->n/2)] = true;

            // Update counts
            match_count[u%(g->n/2)]++;
            match_count[mateI[u]% (g->n/2)]++;

            // Increment the total count
            total_matched_vertices += 2; // Two vertices are matched in each successful match
        }
    }

    // Now, you can print or analyze the match_count vector
    for (int r = 0; r < g->n; ++r) {
        if (match_count[r] > 1) {
            printf("Error: Vertex %d has more than 1 match.\n", r);
        } else if (match_count[r] == 1) {
            //printf("Vertex %d has 1 match.\n", r);
        } else {
            //printf("Vertex %d has 0 matches.\n", r);
        }
    }

    // Print the total number of matched vertices
    printf("Total matched vertices: %d\n", total_matched_vertices);
    printf("Total matched edges: %d\n", total_matched_vertices/2);
    *initial_match_count = total_matched_vertices/2;
    
	//free_graph(g);
	free(g);
	return 0;
}





// Parallel Disjoint BFS with tree grafting

int* MS_BFS_Graft (graph* G, int* mateI)
{
    
   
    
	const int NE = G->m; // number of edges
	const int NV = G->n; // numver of vertices in both sides
    const int nrows = G->nrows; // number of vertices in the left side
	int * restrict endVertex = G->endV; // adjacency
	int * restrict vtx_pointer = G->vtx_pointer; // adjacency pointer
	
	int* QF = (int*) malloc(NV * sizeof(int));
    int* QFnext = (int*) malloc(NV * sizeof(int));
	int* restrict flag = (int*) malloc(NV * sizeof(int));
	int* restrict parent = (int*) malloc(NV * sizeof(int));
	int* restrict leaf = (int*) malloc(NV * sizeof(int));
    int* restrict root = (int*) malloc(NV * sizeof(int));
	int* restrict mate = (int*) malloc(NV * sizeof(int));
    int* unmatchedU = (int*) malloc(nrows * sizeof(int));
    int* nextUnmatchedU = (int*) malloc(nrows * sizeof(int));
    
    
    
  
    double time_start = omp_get_wtime();
    
    #define THREAD_BUF_LEN 16384
    int numUnmatchedU = 0;
    
    // identify unmatched and non-isolated vertices from where search will begin
#pragma omp parallel
    {
        int kbuf=0, nbuf[THREAD_BUF_LEN];
#pragma omp for
        for(int u=0; u<nrows; u++)
        {
            if(mateI[u] == -1 && (vtx_pointer[u+1] > vtx_pointer[u]))
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
                root[u] = u;
            }
            else
                root[u] = -1;
                
                
        }
        if(kbuf>0)
        {
            int voff = __sync_fetch_and_add (&numUnmatchedU, kbuf);
            for (int vk = 0; vk < kbuf; ++vk)
                unmatchedU[voff + vk] = nbuf[vk];
        }
        
    }
    
    
#pragma omp parallel for schedule(static)
    for(int i=0; i<nrows; i++)
    {
        parent[i] = -1;
        leaf[i] = -1;
        mate[i] = mateI[i];
        flag[i] = 0;
    }
    // I separated them out so that root information for row vertices can be set earlier
#pragma omp parallel for schedule(static)
    for(int i=nrows; i<NV; i++)
    {
        parent[i] = -1;
        leaf[i] = -1;
        mate[i] = mateI[i];
        root[i] = -1;
        flag[i] = 0;
    }

    // prepare frontier for the first iteration.
#pragma omp parallel for schedule(static)
	for(int i=0; i<numUnmatchedU; i++) // &&
	{
		int u  = unmatchedU[i];
        QF[i] = u;
        unmatchedU[i] = u;
	}
    
    
    int QFsize = numUnmatchedU;
    int QRsize = NV-nrows;
    int total_aug_path_len=0, total_aug_path_count=0;
    int edgeVisited = 0;
    int eFwdAll=0, eRevAll = 0, eGraftAll=0;
    double timeFwdAll=0, timeRevAll=0, timeGraftAll = 0, timeAugmentAll=0, timeStatAll=0;

    printf("\n************* Starting MS-BFS-Graft Algorithm  *************\n");
    printf("Initial number of non-isolated row vertices = %ld\n\n", numUnmatchedU);

    printf("====================Phase by phase statistics===============================\n");
    printf(" Phase   Initial-unmatched  Matched-in-this-phase    Max-Level     Time (sec)\n");
    printf("============================================================================\n");

    
	int iteration = 1;
	int matched = 1;
	while(matched)
	{
        double time_phase = omp_get_wtime();
        double timeFwd=0, timeRev=0;
        int  phaseEdgeVisited = 0;
        int curLayer = 0;
        int QFnextSize = 0;
        int eFwd=0, eRev = 0;
        int eFwdFrontier = 0;
        double tsLayer;
        
        // Step 1: BFS
#pragma omp parallel
        {
            int kbuf, nbuf[THREAD_BUF_LEN]; // temporary thread private Queue buffer
            while(QFsize > 0)
            {
                bool isTopdownBFS = true;
                double alpha=5;
                if(QFsize > QRsize/alpha)
                    isTopdownBFS=false;

#pragma omp single nowait
                {
                    tsLayer = omp_get_wtime();
                }
                
                
                
                kbuf=0;
                #pragma omp barrier
                //isTopdownBFS=false;
                if(isTopdownBFS) // top-down BFS
                {
                    int teFwd = 0;
#pragma omp for
                    for(int i=0; i<QFsize; i++)
                    {
                        int u = QF[i]; // fairness in Queue acess does not help... tested
                        int curRoot = root[u]; // for unmatched U root[u] = u;
                        if( leaf[curRoot] == -1) // without this test this algorithm is still correct, but continues expanding a dead tree
                        {
                            int j;
                            //for(int j=vtx_pointer[u+1]-1; j>=vtx_pointer[u]; j--)
                            for(j=vtx_pointer[u]; j<vtx_pointer[u+1]; j++)
                            {
                                int v = endVertex[j]; // fairness in accessing adjacenty is not helping. why??: no matter how we access, every neighbor will be in the same tree (in serial case). Hence it does not change #iteration. In DFS this may discover shorter augmenting path, but in BFS it does not help.
                                
                                if(flag[v]==0) // avoids unnecessary __sync_fetch_and_or
                                {
                                    if( __sync_fetch_and_or(&flag[v], 1) == 0 )
                                    {
                                        root[v] = curRoot;
                                        parent[v] = u;
                                        
                                        if(mate[v] == -1)
                                        {
                                            leaf[curRoot] = v; //benign race
                                            break;
                                        }
                                        else
                                        {
                                            int next_u = mate[v];
                                            root[next_u] = curRoot;
                                            
                                            if (kbuf < THREAD_BUF_LEN)
                                            {
                                                nbuf[kbuf++] = next_u;
                                            }
                                            else
                                            {
                                                int voff = __sync_fetch_and_add (&QFnextSize, THREAD_BUF_LEN);
                                                for (int vk = 0; vk < THREAD_BUF_LEN; ++vk)
                                                    QFnext[voff + vk] = nbuf[vk];
                                                nbuf[0] = next_u;
                                                kbuf = 1;
                                            }
                                        }
                                    }
                                }
                            }
                            teFwd += j - vtx_pointer[u];
                            
                        }
                    }
                    __sync_fetch_and_add(&eFwd,teFwd);
                    
                }
                else // bottom up BFS
                {
                    int teRev=0;
#pragma omp for
                    for(int v=nrows; v<NV; v++)
                    {
                        if(flag[v]==0)
                        {
                            int j;
                            for(j=vtx_pointer[v+1]-1; j>=vtx_pointer[v]; j--)  // fairness here is important
                                //for(j=vtx_pointer[v]; j<vtx_pointer[v+1]; j++)
                            {
                                int u= endVertex[j];
                                // u must be in the current layer or current layer+1, both cases are fine
                                // if u in layer+1 we are taking a super step in graph traversal (meaning parent and children can be in Queue)
                                if(root[u]!= -1 && leaf[root[u]] == -1) // u is in an active tree
                                {
                                    // Obtaining a parent in the lowest layer gives the shorter augmenting paths
                                    // But it does not reduce size of frontier at any point
                                    // It requires travesing whole adjaency (without the break below), thus costlier
                                    root[v] = root[u];
                                    parent[v] = u;
                                    flag[v] = 1;
                                    if(mate[v] == -1)
                                    {
                                        leaf[root[v]] = v;  // possible benign race
                                    }
                                    else
                                    {
                                        int next_u = mate[v];
                                        root[next_u] = root[v];
                                        
                                        if (kbuf < THREAD_BUF_LEN)
                                        {
                                            nbuf[kbuf++] = next_u;
                                        }
                                        else
                                        {
                                            int voff = __sync_fetch_and_add (&QFnextSize, THREAD_BUF_LEN);
                                            for (int vk = 0; vk < THREAD_BUF_LEN; ++vk)
                                                QFnext[voff + vk] = nbuf[vk];
                                            nbuf[0] = next_u;
                                            kbuf = 1;
                                        }
                                    }
                                    break;
                                    
                                }
                            }
                            //teRev += j - vtx_pointer[v];
                            teRev += vtx_pointer[v+1] - 1 -j;
                        }
                    }
                    __sync_fetch_and_add(&eRev, teRev);
                    
                }
                if(kbuf>0)
                {
                    int64_t voff = __sync_fetch_and_add (&QFnextSize, kbuf);
                    for (int vk = 0; vk < kbuf; ++vk)
                        QFnext[voff + vk] = nbuf[vk];
                }
                
#pragma omp barrier
#pragma omp single
                {
                    int* t;
                    t = QF;
                    QF = QFnext;
                    QFnext = t;
                    QFsize = QFnextSize;
                    QFnextSize = 0;
                    QRsize = QRsize - QFsize; // underestimate
                    if(isTopdownBFS)
                    {
                        timeFwd += omp_get_wtime() - tsLayer;
                        timeFwdAll += omp_get_wtime() - tsLayer;
                    }
                    else
                    {
                        timeRev += omp_get_wtime() - tsLayer;
                        timeRevAll += omp_get_wtime() - tsLayer;
                    }
                    curLayer ++;
                }
            }
        
        }
        double timeBFS = omp_get_wtime() - time_phase;
        
        // ============================================================
        // ---------- Step2: Augment Matching -------------------------
        // ============================================================
        
        double timeAugment_start = omp_get_wtime();
		int nextUnmatchedSize = 0;
#pragma omp parallel
        {
            int kbuf=0, nbuf[THREAD_BUF_LEN];
            int taug_path_len = 0, taug_path_count = 0;
#pragma omp for
            for(int i=0; i<numUnmatchedU; i++)
            {
                int first_u = unmatchedU[i];
                int last_v = leaf[first_u];
                if(last_v != -1)
                {
                    int v = last_v;
                    taug_path_count++;
                    while(v != - 1)
                    {
                        int u = parent[v];
                        int next_v = mate[u];
                        mate[v] = u;
                        mate[u]=v;
                        v = next_v;
                        taug_path_len += 2;
                    }
                }
                else
                {
                    if (kbuf < THREAD_BUF_LEN)
                    {
                        nbuf[kbuf++] = first_u;
                    }
                    else
                    {
                        int voff = __sync_fetch_and_add (&nextUnmatchedSize, THREAD_BUF_LEN);
                        for (int vk = 0; vk < THREAD_BUF_LEN; ++vk)
                            nextUnmatchedU[voff + vk] = nbuf[vk];
                        nbuf[0] = first_u;
                        kbuf = 1;
                    }
                    //nextUnmatchedU[__sync_fetch_and_add(&nextUnmatchedSize,1)] = first_u;
                }
            }
            if(kbuf>0)
            {
                int voff = __sync_fetch_and_add (&nextUnmatchedSize, kbuf);
                for (int vk = 0; vk < kbuf; ++vk)
                    nextUnmatchedU[voff + vk] = nbuf[vk];
            }
            __sync_fetch_and_add(&total_aug_path_len, taug_path_len);
            __sync_fetch_and_add(&total_aug_path_count, taug_path_count);

            
        }
        
        matched = numUnmatchedU - nextUnmatchedSize; // number of new vertices matched in this phase
       
        int* t = unmatchedU;
		unmatchedU = nextUnmatchedU;
        nextUnmatchedU = t;
        int tNumUnmatchedU = numUnmatchedU;
        numUnmatchedU = nextUnmatchedSize;
		timeAugmentAll  += omp_get_wtime() - timeAugment_start ;
        
        
        
        // ===========================================================================
        // statistics: active, inactive & renewable vertices
        // This statistics are use to decide whether to apply tree grafting mechanism
        // ===========================================================================
        
        double timeStat_start = omp_get_wtime();
        int ActiveVtx = 0, InactiveVtx=0, RenewableVtx = 0;
        int step = 100; // smpling, computing statistics for every 100th vertices
        
#pragma omp parallel
        {
            int tActiveVtx = 0, tInactiveVtx=0, tRenewableVtx = 0; //thread private variables
            
#pragma omp for
            for(int u=0; u<nrows; u+=step)
            {
                if(root[u]!=-1 && leaf[root[u]]==-1)
                    tActiveVtx++;
            }
            
#pragma omp for
            for(int v=nrows; v<NV; v+=step)
            {
                if(root[v]==-1)
                    tInactiveVtx ++;
                else if(leaf[root[v]]!=-1)
                    tRenewableVtx ++;
            }
            
            __sync_fetch_and_add(&ActiveVtx,tActiveVtx);
            __sync_fetch_and_add(&RenewableVtx, tRenewableVtx);
            __sync_fetch_and_add(&InactiveVtx,tInactiveVtx);
        }

        //cout << ActiveVtx << " " << RenewableVtx << " " << InactiveVtx << endl;
        timeStatAll+= omp_get_wtime() - timeStat_start;
        
        
        
        // ===========================================================================
        // Step3: reconstruct frontier for the next iteration
        // Either use tree-grafting or build from scratch
        // ===========================================================================
        
        // decide whther tree-grafting is beneficial or not
        bool isGrafting = true;
        double alpha1=5;
        if(ActiveVtx < RenewableVtx/alpha1)
            isGrafting=false;
        
        isGrafting = true;
        double timeGraft_start = omp_get_wtime(); // store the time to reconstruct frontiers
        int eGraft = 0;
        
        if(isGrafting) // use dead-alive scheme
        {
            
            QFsize = 0;
            QRsize = 0;
#pragma omp parallel
            {
                int teGraft = 0;
                int nbuf[THREAD_BUF_LEN];
                int kbuf = 0;
#pragma omp for
                for(int v = nrows; v < NV; v++)
                {
                    //if( root[v]==-1 || leaf[root[v]]!=-1) // consider both dead and unvisited vertices, enble for testing
                    if( root[v]!=-1 && leaf[root[v]]!=-1) // we will not consider unvisited vertices because they can not be part of frontier
                    {
                        // remove v from the dead tree
                        flag[v] = 0;
                        root[v] = -1;
                        
                        // to obtain best result we look for parents in the adjacenty from high to low indices (given that in forward BFS we traverse adjacenty from low to high indices)
                        int j;
                        for(j=vtx_pointer[v+1]-1; j>=vtx_pointer[v]; j--)
                        //for(j=vtx_pointer[v]; j<vtx_pointer[v+1]; j++)
                        {
                            int u = endVertex[j];
                            if(root[u]!= -1 && leaf[root[u]] == -1) // u in an active tree (no augmenting path found in the latest BFS)
                            {
                                if(mate[v] == -1)
                                {
                                    //  leaf[root[v]] = v;  // we can not allow this because then we will try to destroy a tree that has just found an augmenting path (but the augmentation is yet to perform)!
                                    QF[__sync_fetch_and_add(&QFsize,1)] = u; // insert into the frontier again so that we can discover v; we can insert next_u more than once?? No probem: bottom up does not use QF, top down will explore adjacency only once
                                }
                                else
                                {
                                    root[v] = root[u];
                                    parent[v] = u;
                                    flag[v] = 1;
                                    int next_u = mate[v];
                                    root[next_u] = root[v];
                                    //QF[__sync_fetch_and_add(&QFsize,1)] = next_u; // slow version, shared Queue
                                    if (kbuf < THREAD_BUF_LEN)
                                    {
                                        nbuf[kbuf++] = next_u;
                                    }
                                    else
                                    {
                                        int voff = __sync_fetch_and_add (&QFsize, THREAD_BUF_LEN);
                                        for (int vk = 0; vk < THREAD_BUF_LEN; ++vk)
                                            QF[voff + vk] = nbuf[vk];
                                        nbuf[0] = next_u;
                                        kbuf = 1;
                                    }
                                }
                                break;
                            }
                        }
                        //teGraft += j - vtx_pointer[v];
                        teGraft += vtx_pointer[v+1]-1 - j;
                    }
                    
                }
                
                if(kbuf>0)
                {
                    int voff = __sync_fetch_and_add (&QFsize, kbuf);
                    for (int vk = 0; vk < kbuf; ++vk)
                        QF[voff + vk] = nbuf[vk];
                }
                __sync_fetch_and_add(&eGraft, teGraft);
                
            }
            
            QRsize = NV- nrows - QFsize; // speculative reverse frontier size
            
        }
        else // constructs active trees from scratch
        {
#pragma omp parallel for schedule(static) default(shared)
            for(int v = nrows; v < NV; v++)
            {
                if( root[v]!=-1)
                {
                    flag[v]=0;
                    root[v] = -1;
                    //parent[v] = -1;
                }
            }
            
            
#pragma omp parallel for schedule(static) default(shared)
            for(int v = 0; v < nrows; v++)
            {
                if( root[v]!=-1 && leaf[root[v]]==-1)
                {
                    //flag[v]=0;
                    root[v] = -1; // we need this, otherwise in reverse BFS an Y vertex might attached to a zoombie X vertex
                    //parent[v] = -1;
                }
            }
            
            
            QFsize = numUnmatchedU;
#pragma omp parallel for default(shared)
            for(int i=0; i<QFsize; i++)
            {
                int u = unmatchedU[i];
                QF[i] = u;
                root[u] = u;
            }
            
            QRsize = NV-nrows;
        }
        
        
        
        eRevAll += eRev;
        eFwdAll += eFwd;
        eGraftAll += eGraft;
        phaseEdgeVisited += eRev + eFwd + eGraft;
        edgeVisited += eRev + eFwd + eGraft;
        
        double timeGraft = omp_get_wtime() - timeGraft_start;
        timeGraftAll += timeGraft;
        

        int count = 0;
        for(int i=0; i<G->nrows; i++)
        {
            if(mate[i]==-1) count ++;
        }
		//printf("[%ld]. L=%ld Umnrows=%ld matched= %ld QF=%ld QR=%ld tepsF= %lf tepsR= %lf tFT=%lf eFT=%ld, time=%lf\n",iteration,curLayer,tNumUnmatchedU,matched, QFsize, QRsize, timeFwd, timeRev, timeGraft, eGraft, omp_get_wtime() - time_phase);
        printf("%4ld    %12ld %20ld %18ld %12.2lf\n",iteration, tNumUnmatchedU,matched, curLayer, omp_get_wtime() - time_phase);

        
    iteration++;

        
        
	}
    
	double totalTime = omp_get_wtime() - time_start;
    
  
    //printf("%lf %lf %lf %lf : %lf %ld %lf %lf %ld %lf %lf %lf %lf %ld %lf\n", edgeVisited/1000000.0, totalTime, edgeVisited/(1000000*totalTime), (double)total_aug_path_len/total_aug_path_count, timeFwdAll, eFwdAll, eFwdAll/(1000000*timeFwdAll), timeRevAll, eRevAll, eRevAll/(1000000*timeRevAll),
           //timeAugmentAll, timeStatAll, timeGraft, eGraft, eGraft/(1000000*timeGraft));
    //printf("%ld %lf %lf %ld %lf %lf %lf %lf %lf %lf %lf %ld %ld %ld\n", iteration, (double)100.0*(nrows-numUnmatchedU)/nrows, (double)total_aug_path_len/total_aug_path_count, edgeVisited, edgeVisited/(totalTime*1000000), totalTime, timeFwdAll, timeRevAll, timeAugmentAll, timeStatAll, timeGraftAll,  eFwdAll, eRevAll, eGraft);
	

    // numUnmatchedU contains only non-isolated unmatched vertices
    // compute actual matching cardinality
    int matched_rows = 0;
#pragma omp parallel
    {
        int tmatched = 0; //thread private variables
#pragma omp for
        for(int u=0; u<nrows; u++)
            if(mate[u]!=-1) tmatched++;
        __sync_fetch_and_add(&matched_rows,tmatched);
    }
    
    int isolated_rows = nrows - matched_rows - numUnmatchedU;
    
    printf("============================================================================\n\n");
    printf("========= Overall Statistics ===============\n");
    printf("Number of  Iterations           = %ld \n", iteration);
    printf("Avg. Length of Augmenting Paths = %.2lf \n", (double)total_aug_path_len/total_aug_path_count);
    printf("Total time                      = %.2lf sec\n", totalTime);
    printf("Maximum matching cardinality    = %ld (%.2lf%%)\n", matched_rows*2, (double)100.0*(matched_rows*2)/G->n);
    printf("Matched Rows cardinality        = %ld (%.2lf%%)\n", matched_rows, (double)100.0*(matched_rows)/nrows);
    printf("Isolated Rows                   = %ld (%.2lf%%)\n", isolated_rows, (double)100.0*(isolated_rows)/nrows);
    
    printf("===========================================\n");
    
    printf("============= Time Breakdown ==============\n");
    printf("Total time    = %lf seconds\n", totalTime);
    printf("Top-Down      = %.2lf%%\n", 100*timeFwdAll/totalTime);
    printf("Bottom-UP     = %.2lf%%\n", 100*timeRevAll/totalTime);
    printf("Augmentation  = %.2lf%%\n", 100*timeAugmentAll/totalTime);
    printf("Tree-Grafting = %.2lf%%\n", 100*timeGraftAll/totalTime);
    printf("Statistics    = %.2lf%%\n", 100*timeStatAll/totalTime);
    printf("==========================================\n");
    
    printf("================= Edge Visited ===========\n");
    printf("Total         = %.2lf M\n", edgeVisited/1000000.0);
    printf("Top-Down      = %.2lf M\n", eFwdAll/1000000.0);
    printf("Bottom-UP     = %.2lf M\n", eRevAll/1000000.0);
    printf("Tree-Grafting = %.2lf M\n", eGraftAll/1000000.0);
    printf("===========================================\n");
    //for(int i=0; i<NV; i++)
    //{
    //    mateI[i] = mate[i];
    //}
    
    
	free(flag);
	free(QF);
    free(QFnext);
	free(parent);
	free(leaf);
	free(root);
    free(unmatchedU);
    free(nextUnmatchedU);
    
    return(mate);
}
