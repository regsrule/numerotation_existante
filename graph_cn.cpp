/******************** Includes - Defines ****************/
#include "graph_cn.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>
#include"time_manager.h"
#include"function.h"
#include <x86intrin.h>


Pi_t* pi_allocation(int N)
{
	Pi_t* pi = (Pi_t*)malloc(N*sizeof(int));
	pi->nb_node = N;

	return pi;
}

/***** Memory allocation - graph initialization *****/
Graph_t* graph_allocation(int N, int oriented)
{

	Graph_t* G = (Graph_t*)malloc(sizeof(Graph_t));
	int i;
	G->oriented = oriented;
	G->nb_node  = N;
	G->Nodes    = (Node_t*)malloc(N*sizeof(Node_t));
    
    for (i = 0; i < N; i++)
	{
		G->Nodes[i].deg = 0;
		G->Nodes[i].neighboors = (int*) malloc(sizeof(int));
        }	

    return G;
}

/****** Read graph from txt file edges are in the form
 *	src dst (in each line of the file)            *******/	

void Read_txt(char* filename, Graph_t* G)
{
    
    FILE *fid;

    int src, dst;
    int temp_deg;

    	fid = fopen(filename, "r");
   	if (fid == NULL){printf("Error opening the file\n");}

	while (!feof(fid))
	{
		if (fscanf(fid,"%d\t%d\n", &src, &dst))
		{
		//	printf("%d\t%d\n",src,src);
			G->Nodes[src].deg++;
			temp_deg = G->Nodes[src].deg;
			G->Nodes[src].neighboors = (int*) realloc(G->Nodes[src].neighboors, temp_deg * sizeof(int));
			G->Nodes[src].neighboors[temp_deg - 1] = dst;
		        if(!G->oriented)
			{
				G->Nodes[dst].deg++;
				temp_deg = G->Nodes[dst].deg;
				G->Nodes[dst].neighboors = (int*) realloc(G->Nodes[dst].neighboors, temp_deg * sizeof(int));
				G->Nodes[dst].neighboors[temp_deg - 1] = src;
			}	
				
		}
	}
	//printf("End of connections insertion!\n");
	fclose(fid);
}
