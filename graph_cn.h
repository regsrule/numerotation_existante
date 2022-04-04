#ifndef GRAPH_CN
#define GRAPH_CN

#include <omp.h>
/***** Struct used for Nodes data *****/

typedef struct
{
	int *neighboors;
	int deg;
}Node_t;

typedef struct
{
	int oriented; //0 or 1 for yes or no
	int nb_node;
	Node_t* Nodes;
}Graph_t;

typedef struct
{
	int nb_node;
	int* at;
}Pi_t;

void Read_txt(char* filename, Graph_t* G);
Graph_t* graph_allocation(int N, int oriented);
#endif
