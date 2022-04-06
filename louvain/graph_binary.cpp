// File: graph_binary.cpp
// -- graph handling source
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details
//-----------------------------------------------------------------------------
//Modified by T. Messi Nguele
//For cn-order heuristic, a community-aware ordering for efficient graph analysis
//Jully 2017
//-----------------------------------------------------------------------------

#include <sys/mman.h>
#include <fstream>
#include "graph_binary.h"
#include "math.h"

Graph::Graph() {
  nb_nodes     = 0;
  nb_links     = 0;
  total_weight = 0;
}

Graph::Graph(char *filename, char *filename_w, int type) {
  ifstream finput;
  finput.open(filename,fstream::in | fstream::binary);

  // Read number of nodes on 4 bytes
  finput.read((char *)&nb_nodes, 4);
  assert(finput.rdstate() == ios::goodbit);

  // Read cumulative degree sequence: 8 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
  degrees.resize(nb_nodes);
  finput.read((char *)&degrees[0], nb_nodes*8);

  // Read links: 4 bytes for each link (each link is counted twice)
  nb_links=degrees[nb_nodes-1];
  links.resize(nb_links);
  finput.read((char *)(&links[0]), (long)nb_links*8);  

  // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
  weights.resize(0);
  total_weight=0;
  if (type==WEIGHTED) {
    ifstream finput_w;
    finput_w.open(filename_w,fstream::in | fstream::binary);
    weights.resize(nb_links);
    finput_w.read((char *)&weights[0], (long)nb_links*4);  
  }    

  // Compute total weight
  for (unsigned int i=0 ; i<nb_nodes ; i++) {
    total_weight += (double)weighted_degree(i);
  }
}

Graph::Graph(int n, int m, double t, int *d, int *l, float *w) {
/*  nb_nodes     = n;
  nb_links     = m;
  total_weight = t;
  degrees      = d;
  links        = l;
  weights      = w;*/
}


void
Graph::display() {
/*  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator > p = neighbors(node);
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      if (node<=*(p.first+i)) {
	if (weights.size()!=0)
	  cout << node << " " << *(p.first+i) << " " << *(p.second+i) << endl;
	else
	  cout << node << " " << *(p.first+i) << endl;
      }
    }   
  }*/
  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator > p = neighbors(node);
    cout << node << ":" ;
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      if (true) {
	if (weights.size()!=0)
	  cout << " (" << *(p.first+i) << " " << *(p.second+i) << ")";
	else
	  cout << " " << *(p.first+i);
      }
    }
    cout << endl;
  }
}

void
Graph::display_reverse() {
  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator > p = neighbors(node);
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      if (node>*(p.first+i)) {
	if (weights.size()!=0)
	  cout << *(p.first+i) << " " << node << " " << *(p.second+i) << endl;
	else
	  cout << *(p.first+i) << " " << node << endl;
      }
    }   
  }
}


bool
Graph::check_symmetry() {
  int error=0;
  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator > p = neighbors(node);
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      unsigned int neigh = *(p.first+i);
      float weight = *(p.second+i);
      
      pair<vector<unsigned int>::iterator, vector<float>::iterator > p_neigh = neighbors(neigh);
      for (unsigned int j=0 ; j<nb_neighbors(neigh) ; j++) {
	unsigned int neigh_neigh = *(p_neigh.first+j);
	float neigh_weight = *(p_neigh.second+j);

	if (node==neigh_neigh && weight!=neigh_weight) {
	  cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
	  if (error++==10)
	    exit(0);
	}
      }
    }
  }
  return (error==0);
}


void
Graph::display_binary(char *outfile) {
  ofstream foutput;
  foutput.open(outfile ,fstream::out | fstream::binary);

  foutput.write((char *)(&nb_nodes),4);
  foutput.write((char *)(&degrees[0]),4*nb_nodes);
  foutput.write((char *)(&links[0]),8*nb_links);
}

void
Graph::adjacency_graph2binary_graph(vector<vector<int>> adj_g, char* filename_w, int type) {//added by Messi

  vector<vector<pair<int,float> > > nodes;

  nb_nodes=0;
  int adj_size = adj_g.size(), i = 0;

  while(i < adj_size){
    unsigned int src, dest;
    double weight=1.;
    src = i;
    int nb_neig = adj_g[src].size();
    for(int j = 0; j< nb_neig; j++){
      dest = adj_g[src][j];
      if (nodes.size()<=max(src,dest)+1) {
        nodes.resize(max(src,dest)+1);
	nb_nodes = max(src,dest) + 1;
      }
      
      //if(nodes[src].size()==0)
	     // nb_nodes++;
      nodes[src].push_back(make_pair(dest,weight));
      //We consider that the graph is oriented, otherwise this condition will be used
     /* if (src!=dest)
      {
	if(nodes[dest].size()==0)
	      nb_nodes++;
        nodes[dest].push_back(make_pair(src,weight));
  	nb_links++;
      }*/

      nb_links++;
    }
    i++;
  }

  cout<<"Nb_nodes = "<<nb_nodes<<"\tNb_edges = "<<nb_links<<endl;
  // setting links: each link is counted twice
  // and setting cumulative degrees
  int k = 0;
  long tot=0;
  links.resize(nb_links);
  degrees.resize(nb_nodes);
  for (unsigned int i=0 ; i<nb_nodes ; i++) {
    tot+=(long)nodes[i].size();
    degrees[i] = tot;
    for (unsigned int j=0 ; j<nodes[i].size() ; j++) {
      links[k] =  nodes[i][j].first;
      k++;
    }
  }

  // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
  weights.resize(0);
  total_weight=0;

  // Compute total weight
  for (unsigned int i=0 ; i<nb_nodes ; i++) {
    total_weight += (double)weighted_degree(i);
  }
}

void
Graph::txt_file2binary_graph(char *filename, char* filename_w, int type) {//added by Messi
  ifstream finput;
  finput.open(filename,fstream::in);

  vector<vector<pair<int,float> > > nodes;

  nb_nodes=0;

  while (!finput.eof()) {
    unsigned int src, dest;
    double weight=1.;

      finput >> src >> dest;
    if (finput) {
      if (nodes.size()<=max(src,dest)+1) {
        nodes.resize(max(src,dest)+1);
	nb_nodes = max(src,dest) + 1;
      }
      
      //if(nodes[src].size()==0)
	     // nb_nodes++;
      nodes[src].push_back(make_pair(dest,weight));
      //We consider that the graph is oriented, otherwise this condition will be used
/*if (src!=dest)
      {
	if(nodes[dest].size()==0)
	      nb_nodes++;
        nodes[dest].push_back(make_pair(src,weight));
  	nb_links++;
      }*/

      nb_links++;
    }
  }

  finput.close();
  cout<<"Nb_nodes = "<<nb_nodes<<"\tNb_edges = "<<nb_links<<endl;
  // setting links: each link is counted twice
  // and setting cumulative degrees
  int k = 0;
  long tot=0;
  links.resize(nb_links);
  degrees.resize(nb_nodes);
  for (unsigned int i=0 ; i<nb_nodes ; i++) {
    tot+=(long)nodes[i].size();
    degrees[i] = tot;
    for (unsigned int j=0 ; j<nodes[i].size() ; j++) {
      links[k] =  nodes[i][j].first;
      k++;
    }
  }

  // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
  weights.resize(0);
  total_weight=0;

  // Compute total weight
  for (unsigned int i=0 ; i<nb_nodes ; i++) {
    total_weight += (double)weighted_degree(i);
  }
}
