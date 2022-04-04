#ifndef DEF_FUNCTION
#define DEF_FUNCTION

#define len_string 100 
#define max_char_read 100
#include <stdbool.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include "louvain/graph_binary.h"
#include "louvain/community.h"
#include "rabbit/reorder.h"
#include "gorder/Graph_gorder.h"
#include "gorder/Util.h"

typedef std::vector<std::tuple<int, int, std::vector<int>> > comm_class_t; //community class
typedef std::vector<std::vector<int>> com_h_t; //community hierarchy
typedef std::vector<com_h_t> dendrodram_t; //community hierarchy

void display_time(const char *str);
void print_re_ordered_graph(vector<vector<int>> ReOrderedGraph, const char* filename);
vector <vector<int> > store_neighboors_parallel(vector<vector <int>> adj_graph, vector<int> pi);
vector <vector<int> > store_neighboors_sort(vector<vector <int>> adj_graph, vector<int> pi);
vector <vector<int> > store_neighboors_tiny_sort(vector<vector <int>> adj_graph, vector<int> pi);
vector <vector<int> > store_neighboors(vector<vector <int>> adj_graph, vector<int> pi);
void order_per_community(vector<int> comm_cl, int off_set, int class_size, int comm_min_size, vector<vector <int> >* adj_graph, vector<int> *pi, int W);
//vector<int>gorder_community(Graph_gorder comm_graph, int W);
vector<int>gorder_community(Graph_gorder comm_graph, int W, vector<int> order);
Graph_gorder comm2gorder_graph(vector<vector<int> > adj_graph, int nb_edge);
std::vector<int> sort_communities(std::vector<std::vector<int> > affinity);
comm_class_t classify_and_find_offsets(dendrodram_t tab_comm, int cache_size, vector<vector <int>> adj_graph);
void print_vector(vector<int> vec);
void print_adjacency_graph(vector<vector<int>> com_cl);
void print_communities_class(comm_class_t com_cl);
dendrodram_t detect_communities_louvain_from_graph(vector<vector<int>> adj_g, int type, int nb_pass, float precision, int min_improvement);
dendrodram_t detect_communities_louvain(char* filename, int type, int nb_pass, float precision, int min_improvement);
dendrodram_t detect_communities_rabbit(char* filename); 
char *trim_space(char *str);
char** str_split(char* str, char* car, char*reslt[]);
void print_string(char* elt);
void print_int(int*elt);
#endif
