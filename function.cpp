#include"function.h"

void display_time(const char *str) {
  time_t rawtime;
  time ( &rawtime );
  cerr << str << ": " << ctime (&rawtime);
}
	
	
	
	
vector <vector<int> > store_neighboors_parallel(vector<vector <int>> adj_graph, vector<int> pi)
{
	int nb_node = adj_graph.size(),i;

	vector<vector<int> > graph (adj_graph.size());
  	#pragma omp parallel for schedule(dynamic,1)
	for(i = 0; i < nb_node; i++)
	{
		int nb_neig = adj_graph[i].size(), begin, end;
		vector<int> neighboors(nb_neig);
		for(int j=0; j < nb_neig; j++)
		{
			int neig = adj_graph[i][j];
			neighboors[j] = pi[neig];
		}
		adj_graph[i].clear();
		graph[pi[i]] = neighboors;
	}

	#pragma omp barrier
	return graph;
}

vector <vector<int> > store_neighboors_sort(vector<vector <int>> adj_graph, vector<int> pi)
{
	int nb_node = adj_graph.size(),i;

	vector<vector<int> > graph (adj_graph.size());
	for(i = 0; i < nb_node; i++)
	{
		int nb_neig = adj_graph[i].size(), begin, end;
		vector<int> neighboors(nb_neig);
		for(int j=0; j < nb_neig; j++)
		{
			int neig = adj_graph[i][j];
			neighboors[j] = pi[neig];
		}
		adj_graph[i].clear();
		sort(neighboors.begin(), neighboors.end());
		graph[pi[i]] = neighboors;
	}

	return graph;
}

vector <vector<int> > store_neighboors_tiny_sort(vector<vector <int>> adj_graph, vector<int> pi)
{
	int nb_node = adj_graph.size(),i;

	vector<vector<int> > graph (adj_graph.size());
	for(i = 0; i < nb_node; i++)
	{
		int nb_neig = adj_graph[i].size(), begin, end;
		vector<int> neighboors(nb_neig);
		begin = 0;
		end   = nb_neig-1;
		//cout<<"i= "<<i<<endl;
		for(int j=0; j < nb_neig; j++)
		{
			int neig = adj_graph[i][j];
			if(pi[i] > pi[neig])
			{
				neighboors[begin] = pi[neig];
				begin++;
			//cout<<"\tj= "<<j<<" begin ="<<begin<<" end = "<<end<<endl;
			}
			else
			{
				neighboors[end] = pi[neig];
				end--;
			//cout<<"\tj= "<<j<<" begin ="<<begin<<" end = "<<end<<endl;
			}
		}
		adj_graph[i].clear();
		graph[pi[i]] = neighboors;
	}

	return graph;
}

vector <vector<int> > store_neighboors(vector<vector <int>> adj_graph, vector<int> pi)
{
	int nb_node = adj_graph.size(),i;

	vector<vector<int> > graph (adj_graph.size());
	for(i = 0; i < nb_node; i++)
	{
		int nb_neig = adj_graph[i].size(), begin, end;
		vector<int> neighboors(nb_neig);
		for(int j=0; j < nb_neig; j++)
		{
			int neig = adj_graph[i][j];
			neighboors[j] = pi[neig];
		}
		adj_graph[i].clear();
		graph[pi[i]] = neighboors;
	}

	return graph;
}

Graph_gorder comm2gorder_graph(vector<vector<int> > adj_graph, int nb_edge)
{
	clock_t start, end;
	string filename;
	Graph_gorder g;
	string name;
	srand(time(0));

	start=clock();
	//cout << "comm_size = "<<adj_graph.size()<<endl;
	g.setGraph(adj_graph, nb_edge);
	g.Transform();
	//cout << "\tcomm_size = "<<adj_graph.size()<<endl;
/*
	cout << name << " building community graph is complete." << endl;
	end=clock();
	cout << "Time Cost: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
*/
	return g;
}

vector<int>gorder_community(Graph_gorder comm_graph, int W, vector<int> order)
{
	clock_t start, end;

	start=clock();
	comm_graph.GorderGreedy(order, W);
	end=clock();

	/*cout << "ReOrdered Time Cost: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
	cout << "Begin Output the Reordered Graph" << endl;
	comm_graph.PrintReOrderedGraph(order);
	cout << endl;
*/
	return order;
}

void order_per_community(vector<int> comm_cl, int off_set, int class_size, int comm_min_size, vector<vector <int> >* adj_graph, vector<int> *pi, int W)
{
	clock_t start, end, nb_edge = 0;
	int comm_size = comm_cl.size();
       	vector <int> result_order;
	vector<vector<int> > comm_graph(comm_cl.size());

	if(class_size <= comm_min_size)
	{
		for(int i = 0; i<comm_size; i++)
		{
			int old_pos = comm_cl[i];   //the initial position of the node
			(*pi)[old_pos] = i + off_set;  //setting the final position of the node
		}
		return;
	}

	vector<int> simple_order(comm_size);
	vector<int> graph_map((*adj_graph).size(),-1);
	for(int i=0; i<comm_size; i++)
	{
		simple_order[i] = comm_cl[i];
		graph_map[comm_cl[i]] = i;
	}

	//preparing adjacency list community used to convert into gorder graph
	int i = 0, r_nodes = 0, not_r_nodes = 0;
	while(i<comm_size)
	{
		vector<int> neighs;
		int new_nb_neig, nb_neig, real_node, real_neig;
		new_nb_neig = 0;
		real_node = comm_cl[i];
		nb_neig     = (*adj_graph)[real_node].size();
		neighs.reserve(nb_neig);
		for(int j = 0; j<nb_neig; j++)
		{
			real_neig = (*adj_graph)[real_node][j];
			if(graph_map[real_neig]!=-1)
			{
				neighs.push_back(graph_map[real_neig]);
				nb_edge++;
				new_nb_neig++;
			}
		}

		comm_graph[i] = neighs;
		i++;
	}
	

	Graph_gorder    g = comm2gorder_graph(comm_graph, nb_edge);
	vector<int> order;
	
	order = gorder_community(g, W, order);

	for(i = 0; i<order.size();i++)
	{
		int old_pos = simple_order[order[i]]; //the initial position of the node

		(*pi)[old_pos] = i + off_set;            //setting the final position of the node
	}

	return;
}

vector<int> sort_communities(vector<vector<int> > affinity)
{
	int i, j, k, ind_max, max, nb_comm = affinity.size();
	vector <int> successor(nb_comm);
	vector <bool> taken(nb_comm,false);	
	bool found;

	k        = 0; 
	i        = 0; 
	ind_max  = -1;
	max      = -1;
	found    = false;
	taken[i] = true;

	while(k<nb_comm){
		for(j=0;j<nb_comm;j++){
			if(!taken[j] && affinity[i][j] > max){
				ind_max = j;
				max     =  affinity[i][j];
				found   = true;
			}
		}
		if(found){
			successor[i]   = ind_max;
			taken[ind_max] = true;
			i              = ind_max;
		}
		found   = false;
		ind_max = -1;
		max     = -1;
		k ++;
	}
	return successor;
}

comm_class_t classify_and_find_offsets(dendrodram_t tab_comm, int cache_size, vector<vector <int>> adj_graph)
{
	comm_class_t com_cl;
	int nb_hierarchy = tab_comm.size();
	vector <vector <int> > affinity = tab_comm[nb_hierarchy -1];
	unsigned long int k, nb_comm = affinity.size(), cl_num = 0, offset = 0;
	int h;
	vector <int> comm_list;
	vector <int> successor;

	successor = sort_communities(affinity);

	
	int c = 0, comm = 0;
	while(comm < nb_comm)
	{
		h = nb_hierarchy - 2;
		comm_list = tab_comm[h][c];
		if(h == 0)//We have directly communities without subcommunities
		{
			com_cl.resize(cl_num+1);
			unsigned long int nb_nodes_sc = tab_comm[h][c].size(), class_size = 0;
			unsigned long int base = std::get<2>(com_cl[cl_num]).size();
			std::get<2>(com_cl[cl_num]).resize(base + nb_nodes_sc);
			for(k = 0; k < nb_nodes_sc; k++)
			{
			   int node = tab_comm[h][c][k];
			   class_size += adj_graph[node].size();
			   std::get<2>(com_cl[cl_num])[base + k] = node;
			}
			std::get<0>(com_cl[cl_num]) = offset;
			std::get<1>(com_cl[cl_num]) = class_size;
			offset = offset + std::get<2>(com_cl[cl_num]).size();
			cl_num++;
			class_size = 0;
		}
		else
		{
			h = h-1;
			do{
			   if(h == 0)
			   {
				   unsigned int sc = 0, class_size = 0;
				   vector <int> n_comm_list;
				   vector <int> pn_comm_list;
				   com_cl.resize(cl_num+1);
				   while(sc < comm_list.size())
				   {
					   unsigned long int nb_nodes_sc = tab_comm[h][comm_list[sc]].size();
					   unsigned long int base = std::get<2>(com_cl[cl_num]).size();
					   std::get<2>(com_cl[cl_num]).resize(base + nb_nodes_sc);
					   for(k = 0; k < nb_nodes_sc; k++)
					   {
						   int node = tab_comm[h][comm_list[sc]][k];
						   class_size += adj_graph[node].size();
						   std::get<2>(com_cl[cl_num])[base + k] = node;
						   //if(nb_nodes_sc == 1 || nb_nodes_sc == 2)
						  //	 cout<<"c "<<c<<"\tclass "<<cl_num<<"\tsc = "<<sc<<"\tnb_nodes_sc = "<<nb_nodes_sc<<"\t node = "<<node<<endl;
					   }

					   if((class_size >= cache_size) || (sc == (comm_list.size()-1)))
					   {
						std::get<0>(com_cl[cl_num]) = offset;
						std::get<1>(com_cl[cl_num]) = class_size;
						offset = offset + std::get<2>(com_cl[cl_num]).size();
						cl_num++;
						class_size = 0;
						if(sc < comm_list.size())
						   com_cl.resize(cl_num+1);

					   }
					sc++;
				   }
			   }
			   else
			   {
				   vector <int> n_comm_list;
				   unsigned int sc = 0;
				   while(sc < comm_list.size())
				   {
					   unsigned long int nb_nodes_sc = tab_comm[h][comm_list[sc]].size();
					   unsigned long int base = n_comm_list.size();
					   n_comm_list.resize(base + nb_nodes_sc);
					   for(k = 0; k < nb_nodes_sc; k++)
					   {
						   int node = tab_comm[h][comm_list[sc]][k];
						   n_comm_list[base + k] = node;
					   }
					sc++;
				   }
				   comm_list = n_comm_list;
			   }
			   h = h-1;
			} while(h>=0);
		}
		c = successor[c];
		comm++;
	}

	
	return com_cl;
}

void print_vector(vector<int> vec)
{
	for(int i = 0; i<vec.size(); i++)
	{
		cerr << i <<":"<< vec[i]<<endl;
	}
}
void print_re_ordered_graph(vector<vector<int>> ReOrderedGraph, const char* filename){
	ofstream out(filename);
	cout<<"Printing reodered graph in file "<<filename<<endl;
	int vsize = ReOrderedGraph.size();
	for(int u=0; u<vsize; u++){
		for(int j=0; j<ReOrderedGraph[u].size(); j++){
			out << u << '\t' << ReOrderedGraph[u][j] << endl;
		}
	}
	out.close();
}

void print_adjacency_graph(vector<vector<int>> com_cl)
{
	cerr << "Printing graph ..." << endl;
	for(int cl = 0; cl < com_cl.size(); cl++)
	{
		cerr << cl <<"("<<com_cl[cl].size() <<")" << ":" ;
		for(int sc = 0; sc < com_cl[cl].size(); sc++)
			cerr << com_cl[cl][sc] << "\t" ;
		cerr << endl;
	}
}

void print_communities_class(comm_class_t com_cl)
{
	for(int cl = 0; cl < com_cl.size(); cl++)
	{
		cerr << cl <<"("<<std::get<0>(com_cl[cl]) <<","<<std::get<1>(com_cl[cl])<<")" << ":" ;
		for(int sc = 0; sc < std::get<2>(com_cl[cl]).size(); sc++)
			cerr << std::get<2>(com_cl[cl])[sc] << "\t" ;
		cerr << endl;
	}
}

dendrodram_t detect_communities_louvain_from_graph(vector<vector<int>> adj_g, int type, int nb_pass, float precision, int min_improvement) 
{
	char* filename_w = NULL, *filename_part = NULL;
 	Community c(adj_g, type, nb_pass, precision);
	dendrodram_t tab_comm;

	bool improvement=true, verbose = false;
	double mod=c.modularity(), new_mod;
	int level=0, i= 0;
	int display_level = -2;
	time_t time_begin, time_end;
	time(&time_begin);
	
	vector <vector<int> > comm_affinity;
	bool store_adj_g = true;

	do {
		Graph g;
		improvement = c.one_level(min_improvement);
		new_mod = c.modularity();

		if(c.size <=10000 && i>1)
			tab_comm.push_back(c.store_partition_louvain(true));
		else
			tab_comm.push_back(c.store_partition_louvain(false));
				
		comm_affinity = c.comm_affinity;
		g = c.partition2graph_binary();
		c = Community(g, nb_pass, precision);

		mod=new_mod;

		if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
		improvement=true;
		i++;
	} while(improvement);

	time(&time_end);
	cout << new_mod << endl;

	tab_comm.push_back(comm_affinity);
	tab_comm.push_back(adj_g);
	return tab_comm;
}

dendrodram_t detect_communities_louvain(char* filename, int type, int nb_pass, float precision, int min_improvement) 
{
	char* filename_w = NULL, *filename_part = NULL;
 	Community c(filename, type, nb_pass, precision);
	vector <vector<int> > adj_g;
	dendrodram_t tab_comm;

	bool improvement=true, verbose = false;
	double mod=c.modularity(), new_mod;
	int level=0, i= 0;
	int display_level = -2;
	time_t time_begin, time_end;
	time(&time_begin);
	
	vector <vector<int> > comm_affinity;
	bool store_adj_g = true;

	adj_g = c.adjacency_graph;
	do {
		Graph g;
		improvement = c.one_level(min_improvement);
		new_mod = c.modularity();

		if(c.size <=10000 && i>1)
		{
			tab_comm.push_back(c.store_partition_louvain(true));
			cout<<"storing affinities ..."<<endl;
		}
		else
			tab_comm.push_back(c.store_partition_louvain(false));
				
		comm_affinity = c.comm_affinity;
		g = c.partition2graph_binary();
		c = Community(g, nb_pass, precision);

		mod=new_mod;

		if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
		improvement=true;
		i++;
	} while(improvement);

	time(&time_end);
	cout << new_mod << endl;

	tab_comm.push_back(comm_affinity);
	tab_comm.push_back(adj_g);
	return tab_comm;
}

dendrodram_t detect_communities_rabbit(char* filename) 
{
  using boost::adaptors::transformed;
	dendrodram_t tab_comm;
	std::vector<std::vector<int> > comm_affinity;
  const std::string graphpath = filename;

  std::cerr << "Number of threads: " << omp_get_max_threads() << std::endl;

  std::cerr << "Reading an edge-list file: " << graphpath << std::endl;
  auto       adj = read_graph(graphpath);
  const auto m   =
      boost::accumulate(adj | transformed([](auto& es) {return es.size();}),
                        static_cast<size_t>(0));
  std::cerr << "Number of vertices: " << adj.size() << std::endl;
  std::cerr << "Number of edges: "    << m          << std::endl;

    tab_comm = store_partition_rabbit(std::move(adj), true);

    return tab_comm;
}

void print_int(int*elt)
{
	printf("%d",*elt);
}

void print_string(char* elt)
{
	printf("%s",elt);
}

char *trim_space(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

char** str_split(char* str, char* car, char*reslt[])
{
	char *str1, *str2, *token, *subtoken;
        char *saveptr1, *saveptr2;
        int j;

	for (j = 1, str1 = str; ; j++, str1 = NULL) {
               token = strtok_r(str1, car, &saveptr1);
               if (token == NULL)
                   break;
		reslt[j-1] = token;
          }
	
	return reslt;
}
