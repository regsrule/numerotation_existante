#include<stdio.h>
#include <limits.h>
#include <sys/time.h>
#include "graph_cn.h"
#include "function.h"
#include <unistd.h>

int else_where_cache_size(int cache_size)
{
	return cache_size/sizeof(int);
}

int first_level_dcache_size()
{
	long l1_cache_size = sysconf(_SC_LEVEL1_DCACHE_SIZE);

	printf("L1 Cache Size is %ld bytes. (%d int)\n", l1_cache_size, l1_cache_size/sizeof(int)); 

	return l1_cache_size / sizeof(int);
}

int second_level_cache_size()
{


	long l2_cache_size = sysconf(_SC_LEVEL2_CACHE_SIZE);
	long l3_cache_size = sysconf(_SC_LEVEL3_CACHE_SIZE);
	long l4_cache_size = sysconf(_SC_LEVEL4_CACHE_SIZE);

	printf("Choosen: L2 Cache Size is %ld bytes. (%d int) \n", l2_cache_size, l2_cache_size/sizeof(int)); 
	printf("L3 Cache Size is %ld bytes.(%d int)\n", l3_cache_size, l3_cache_size/sizeof(int) ); 
	printf("L4 Cache Size is %ld bytes. (%d int)\n", l4_cache_size,l4_cache_size/sizeof(int)); 

	return l2_cache_size / sizeof(int);
}

int last_level_cache_size()
{


	long l2_cache_size = sysconf(_SC_LEVEL2_CACHE_SIZE);
	long l3_cache_size = sysconf(_SC_LEVEL3_CACHE_SIZE);
	long l4_cache_size = sysconf(_SC_LEVEL4_CACHE_SIZE);

	printf("L2 Cache Size is %ld bytes. (%d int) \n", l2_cache_size, l2_cache_size/sizeof(int)); 
	printf("L3 Cache Size is %ld bytes.(%d int)\n", l3_cache_size, l3_cache_size/sizeof(int) ); 
	printf("L4 Cache Size is %ld bytes. (%d int)\n", l4_cache_size,l4_cache_size/sizeof(int)); 

	if(l4_cache_size != 0)
		return l4_cache_size / sizeof(int);
	else
		if(l3_cache_size!=0)
			return l3_cache_size / sizeof(int);
		else
			return l2_cache_size / sizeof(int);
}

vector <vector <int> > g_order(char* graph_file, int W)
{

  	const double tstart1 = rabbit_order::now_sec();
	ios::sync_with_stdio(false);
	int i;
	clock_t start, end;
	string filename = graph_file;

	srand(time(0));

	Graph_gorder g;
	string name;
  	const double tstart4 = rabbit_order::now_sec();
	name=extractFilename(filename.c_str());
	g.setFilename(name);

	start=clock();
	g.readGraph(filename);
	g.Transform();
        cout<<"Reading graph: " << rabbit_order::now_sec() - tstart4 << std::endl;

	start=clock();
	vector<int> order;
  	const double tstart3 = rabbit_order::now_sec();
	g.GorderGreedy(order, W);
        cout<<"Ordering: " << rabbit_order::now_sec() - tstart3 << std::endl;

  	const double tstart2 = rabbit_order::now_sec();
	vector <vector <int> > adj_graph = g.SetReOrderedGraph(order);
        cout<<"Storing new graph: " << rabbit_order::now_sec() - tstart2 << std::endl;

	cout<<"Total time: " << rabbit_order::now_sec() - tstart1 << std::endl;
	return adj_graph;
}

vector <vector <int> > r_order(char* graph_file)
{
	using boost::adaptors::transformed;
	vector <vector <int> > adj_graph;

  	const double tstart1 = rabbit_order::now_sec();
	std::cerr << "Number of threads: " << omp_get_max_threads() << std::endl;

	std::cerr << "Reading an edge-list file: " << graph_file << std::endl;
	auto       adj = read_graph(graph_file);
	const auto m   =
	boost::accumulate(adj | transformed([](auto& es) {return es.size();}),
			static_cast<size_t>(0));
	std::cerr << "Number of vertices: " << adj.size() << std::endl;
	std::cerr << "Number of edges: "    << m          << std::endl;

  	const double tstart2 = rabbit_order::now_sec();
	adj_graph = reorder_set_graph(std::move(adj));
        cout<<"Storing graph: " << rabbit_order::now_sec() - tstart2 << std::endl;

	cout<<"Total time: " << rabbit_order::now_sec() - tstart1 << std::endl;
	return adj_graph;
}

vector <vector <int> > cn_order_rabbit(char* graph_file, int cache_size, int comm_min_size, int W)
{
  	const double tstart1 = rabbit_order::now_sec();
	dendrodram_t tab_comm_r = detect_communities_rabbit(graph_file);
        cout<<"Detecting comm and affinities: " << rabbit_order::now_sec() - tstart1 << std::endl;

	vector <vector <int> > adj_graph_r = tab_comm_r[tab_comm_r.size()-1];
	vector <vector <int> > result_graph;
	vector<int> pi(adj_graph_r.size(), -1);

	tab_comm_r.resize(tab_comm_r.size()-1);
  	const double tstart2 = rabbit_order::now_sec();
	comm_class_t  comm_cl_r  = classify_and_find_offsets(tab_comm_r, cache_size, adj_graph_r);
        cout<<"Classifying: " << rabbit_order::now_sec() - tstart2 << std::endl;


  	const double tstart3 = rabbit_order::now_sec();
  	#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i<comm_cl_r.size(); i++)
	{
		order_per_community(std::get<2>(comm_cl_r[i]), std::get<0>(comm_cl_r[i]), std::get<1>(comm_cl_r[i]),comm_min_size, &adj_graph_r, &pi, W);
		std::get<2>(comm_cl_r[i]).clear();
	}
	#pragma omp barrier
        cout<<"Ordering per community: " << rabbit_order::now_sec() - tstart3 << std::endl;

//	print_vector(pi);
  	const double tstart4 = rabbit_order::now_sec();
	//result_graph = store_neighboors(adj_graph_r, pi);
	result_graph = store_neighboors_sort(adj_graph_r, pi);
        cout<<"Storing neighboors: " << rabbit_order::now_sec() - tstart4 << std::endl;
        
	cout<<"Total time: " << rabbit_order::now_sec() - tstart1 << std::endl;

	return result_graph;
}

vector <vector <int> > cn_order_louvain(char* graph_file, int cache_size, int comm_min_size, int W, int nb_pass, int nb_moves)
{
  	const double tstart1 = rabbit_order::now_sec();
	dendrodram_t tab_comm_l = detect_communities_louvain(graph_file, 1, nb_pass, 0.0001, nb_moves);
        cout<<"Detecting comm and affinities: " << rabbit_order::now_sec() - tstart1 << std::endl;
	vector <vector <int> > adj_graph_l = tab_comm_l[tab_comm_l.size()-1];
	tab_comm_l.resize(tab_comm_l.size()-1);

  	const double tstart2 = rabbit_order::now_sec();
	comm_class_t  comm_cl_l  = classify_and_find_offsets(tab_comm_l, cache_size, adj_graph_l);
        cout<<"Classifying: " << rabbit_order::now_sec() - tstart2 << std::endl;

	//print_communities_class(comm_cl_l);
	vector<int> pi(adj_graph_l.size(), -1);

//	print_adjacency_graph(adj_graph_l); 
  	const double tstart3 = rabbit_order::now_sec();
  	std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
  	#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i<comm_cl_l.size(); i++)
	{
		//if((i+1)<comm_cl_l.size()&& std::get<0>(comm_cl_l[i])<=2573504 && std::get<0>(comm_cl_l[i+1])>=253504)
//cout<<"comm = "<<i<<"\toffset "<<std::get<0>(comm_cl_l[i])<<"\t size "<<std::get<1>(comm_cl_l[i]) <<"\tnext offset "<<std::get<0>(comm_cl_l[i+1])<<"\tNext size "<<std::get<1>(comm_cl_l[i+1])<<endl;
		order_per_community(std::get<2>(comm_cl_l[i]), std::get<0>(comm_cl_l[i]), std::get<1>(comm_cl_l[i]),comm_min_size, &adj_graph_l, &pi, W);
		std::get<2>(comm_cl_l[i]).clear();
	}
        cout<<"Ordering per community: " << rabbit_order::now_sec() - tstart3 << std::endl;

//	print_vector(pi);
	
	#pragma omp barrier
  	const double tstart4 = rabbit_order::now_sec();
	//vector<vector <int>> result_graph = store_neighboors(adj_graph_l, pi);
	vector<vector <int>> result_graph = store_neighboors_sort(adj_graph_l, pi);
        cout<<"Storing neighboors: " << rabbit_order::now_sec() - tstart4 << std::endl;

	cout<<"Total time: " << rabbit_order::now_sec() - tstart1 << std::endl;
	return result_graph; 
}

vector <vector <int> > numbaco_order(vector<vector<int>> adj_graph, char*graph_file, bool is_file, int nb_pass, int min_moves, int cache_size, int comm_min_size, int W)
{

	/*int comm_min_size =  first_level_dcache_size();
	int cache_size    =  last_level_cache_size(), W = 3;*/

  	const double tstart1 = rabbit_order::now_sec();
	dendrodram_t tab_comm_l;
	if(is_file)
		tab_comm_l = detect_communities_louvain(graph_file, 1, nb_pass, 0.0001, min_moves);
	else
	{
		//print_adjacency_graph(adj_graph); 
		tab_comm_l = detect_communities_louvain_from_graph(adj_graph, 1, nb_pass, 0.0001, min_moves);
	}
		
        cout<<"Detecting comm and affinities: " << rabbit_order::now_sec() - tstart1 << std::endl;
	vector <vector <int> > adj_graph_l = tab_comm_l[tab_comm_l.size()-1];
	tab_comm_l.resize(tab_comm_l.size()-1);
	//print_adjacency_graph(adj_graph_l); 

  	const double tstart2 = rabbit_order::now_sec();
	comm_class_t  comm_cl_l  = classify_and_find_offsets(tab_comm_l, cache_size, adj_graph_l);
        cout<<"Classifying: " << rabbit_order::now_sec() - tstart2 << std::endl;

	//print_communities_class(comm_cl_l);
	vector<int> pi(adj_graph_l.size(), -1);

//	print_adjacency_graph(adj_graph_l); 
  	const double tstart3 = rabbit_order::now_sec();
  	std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
  	#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i<comm_cl_l.size(); i++)
	{
		//order_per_community(std::get<2>(comm_cl_l[i]), std::get<0>(comm_cl_l[i]), std::get<1>(comm_cl_l[i]),comm_min_size, &adj_graph_l, &pi, W);
		order_per_community(std::get<2>(comm_cl_l[i]), std::get<0>(comm_cl_l[i]), comm_min_size, comm_min_size, &adj_graph_l, &pi, W);//to avoid gorder class_size = comm_min_size
		std::get<2>(comm_cl_l[i]).clear();
	}
        cout<<"Ordering per community: " << rabbit_order::now_sec() - tstart3 << std::endl;

	//print_vector(pi);
	
	#pragma omp barrier
  	const double tstart4 = rabbit_order::now_sec();
	//vector<vector <int>> result_graph = store_neighboors(adj_graph_l, pi);
	vector<vector <int>> result_graph = store_neighboors_sort(adj_graph_l, pi);
        cout<<"Storing neighboors: " << rabbit_order::now_sec() - tstart4 << std::endl;

	cout<<"Total time: " << rabbit_order::now_sec() - tstart1 << std::endl;
	//print_adjacency_graph(result_graph); 
	return result_graph; 
}

vector <vector <int> > cn_order_naive(char* graph_file, int nb_pass, int min_moves,int cache_size, int comm_min_size, int W)
{
  	const double tstart1 = rabbit_order::now_sec();
	vector<vector <int>> adj_graph = g_order(graph_file, W);


	vector<vector <int>> result_graph = numbaco_order(std::move(adj_graph),NULL,false, nb_pass, min_moves, cache_size, comm_min_size, W);
        cout<<"Total cn_naive: " << rabbit_order::now_sec() - tstart1 << std::endl;

	return result_graph; 
}

int main(int argc, char* argv[])
{
	int comm_min_size, cache_size;
	int             W = 5;

	if (argc != 3)
	{
			cout<<"error not called properly, \ncall as follow: \n ./cn-order filename [0-5], respectively \n cn-order_l ,\n cn-order_r,\n cn-order_naive,\n numbaco_order,\n rabbit_order,\n g_order \n "<<endl;
		
		return -1;

	}



	int algo = atoi(argv[2]);
	//int nb_pass = atoi(argv[3]), min_moves = argv[4]);
	int nb_pass = 10, min_moves = 0;
	vector <vector <int> > adj_graph; 
	
	//print_adjacency_graph(adj_graph_l); 
	
/*	print_adjacency_graph(adj_graph_r); */
	if(algo<=3){

	/*	comm_min_size =  else_where_cache_size(32768);
		cache_size    =  else_where_cache_size(2097152);*/
		comm_min_size =  first_level_dcache_size();
		cache_size    =  last_level_cache_size();
		//cache_size    =  second_level_cache_size();
		W = 5;
	}
	string filename = extractFilename(argv[1]);
	switch(algo)
	{
		case 0:
		case 1:
			if(algo == 0){	
				adj_graph = cn_order_louvain(argv[1], cache_size, comm_min_size, W, nb_pass, min_moves);
				print_re_ordered_graph(adj_graph, (filename+"_cn_louvain.txt").c_str());
			}
			else{
				adj_graph = cn_order_rabbit(argv[1], cache_size, comm_min_size, W);
				print_re_ordered_graph(adj_graph, (filename+"_cn_rabbit.txt").c_str());
			}
			break;
	
		case 2: adj_graph = cn_order_naive(argv[1], nb_pass, min_moves, cache_size, comm_min_size, W);
			print_re_ordered_graph(adj_graph, (filename+"_cn_naive.txt").c_str());
			break;

		case 3: adj_graph = numbaco_order(adj_graph,argv[1],true, nb_pass, min_moves, cache_size, comm_min_size, W);
			print_re_ordered_graph(adj_graph, (filename+"_numbaco.txt").c_str());
			break;

		case 4: adj_graph = r_order(argv[1]);
			print_re_ordered_graph(adj_graph, (filename+"_rabbit.txt").c_str());
			break;

		case 5: adj_graph = g_order(argv[1],W);
			print_re_ordered_graph(adj_graph, (filename+"_gorder.txt").c_str());
			break;

		default:
			cout<<"./cn-order filename [0-5], respectively cn-order_l, cn-order_r, cn-order_naive, numbaco_order, rabbit_order, g_order"<<endl;
	}

	return 0;
}

