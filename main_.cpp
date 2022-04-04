#include<stdio.h>
#include <limits.h>
#include <sys/time.h>
#include "graph_cn.h"


void cn_order(Graph_t* G, Pi_t* pi, int cache_size)
{
	int i, nb_class;

/*  	Com    	 = detect_communities(G);
	Com_cl 	 = classify_find_offsets(Com, cache_size);

	nb_class = Com_cl.nb_class;
#pragma omp parallel for schedule(dynamic,1)
	for(i=0; i<nb_class; i++)
		compute_community_order(Com_cl[i],G,pi);
	
	store_neighboors(G,pi);*/

	return;
}

int main(int argc, char* argv[])
{
	return 0;
}
