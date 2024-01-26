#ifndef _graphgenBP_h
#define _graphgenBP_h



typedef struct /* the bipartite graph data structure */
{
	int n; // numver of vertices in both sides
	int nrows; // number of vertices in the left side
	int m; // number of edges
	int* vtx_pointer; // an array of size n+1 storing the pointer in endV array
    int* endV; //an array of size m that stores the second vertex of an edge.
	double* weight; // not used in unweighted graph
} graph;



void process_mtx_compressed(char *fname, graph* bGraph);
void process_mtx_compressed(char *fname, graph* bGraph, int **rows, int **cols, int **matching, int*nr_ptr, int*nc_ptr, int*nn_ptr,double *parse_graph_time, double *create_csr_time);
void fast_mtx_read_build(char *fname, graph* bGraph);
bool isEqual(graph* bGraph1, graph* bGraph2);

/** Clean up. */
void free_graph (graph* bGraph);
graph* swap_side(graph* bGraph);

#endif
