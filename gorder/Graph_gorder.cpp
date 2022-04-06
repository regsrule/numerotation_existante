/*
MIT License

Copyright (c) 2016, Hao Wei.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//-----------------------------------------------------------------------------
//Modified by T. Messi Nguele
//For cn-order heuristic, a community-aware ordering for efficient graph analysis
//Jully 2017
//-----------------------------------------------------------------------------

#include "Graph_gorder.h"

string Graph_gorder::getFilename(){
	return name;
}

void Graph_gorder::setFilename(string name){
	this->name.swap(name);
}

Graph_gorder::Graph_gorder() {
	edgenum=vsize=0;
}

Graph_gorder::~Graph_gorder() {
}

void Graph_gorder::clear() {
	vsize = 0;
	edgenum=0;
	name.clear();
	graph.clear();
	outedge.clear();
	inedge.clear();
}


void Graph_gorder::TransformSettedGraph(){//Is it time consuming, if so we should change it
	vector<int> order;
	RCMOrder(order);
	if(order.size()!=vsize){
		cout << "order.size()!=vsize" << endl;
		quit();
	}
	if(graph.size()!=(vsize+1)){
		cout << "graph.size()!=(vsize+1)" << endl;
		quit();
	}

	int count_edge = 0;
	vector<int>().swap(inedge);
	vector< pair<int, int> > edges(edgenum);
	//edges.reserve(edgenum);
	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart, limit=graph[i+1].outstart; j<limit; j++)
		{
			edges[count_edge].first  = order[i];
			edges[count_edge].second = order[outedge[j]];
			count_edge++;
			//edges.push_back(make_pair(order[i], order[outedge[j]]));
		}
	}
	if(edges.size()!=edgenum){
		cout << "edges.size()!=edgenum" << endl;
		quit();
	}

	for(int i=0; i<vsize; i++){
		graph[i].outdegree=graph[i].indegree=0;
	}
	for(int i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}

	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;

	std::sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{

			if(a.second<=b.second)
				return true;
			else
				return false;
		}
	});

	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}
	vector< pair<int, int> >().swap(edges);
	vector<int> inpos(vsize);
	for(int i=0; i<vsize; i++){
		inpos[i]=graph[i].instart;
	}
	inedge.resize(edgenum);
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outstart+graph[u].outdegree; j++){
			inedge[inpos[outedge[j]]]=u;
			inpos[outedge[j]]++;
		}
	}
}

void Graph_gorder::setGraph(vector<vector<int>> adj_graph, int nb_edge){

	char line[40];
	int u, v, _vsize = 0;
	const char* str=NULL;
	
	vsize=0;
	edgenum=0;
	vector< pair<int, int> > edges;
	vector<bool> v_present(adj_graph.size(),false);
	//edges.reserve(100000000);
	edges.reserve(nb_edge);
	int i = 0;
	while(i<adj_graph.size()){
			u= i;
			if(u == 31147 && adj_graph.size() == 31148)
				cout<<"adj_graph["<<i<<"].size() = "<<adj_graph[i].size()<<endl;
		for(int j = 0; j < adj_graph[i].size(); j++){
			v=adj_graph[u][j];
			if(u == 31147 && adj_graph.size() == 31148)
				cout<<"Normally something should be donne"<<endl;
			if(u==v)
			{
				cout<<u<<"("<<adj_graph.size()<<")"<<endl;
				continue;
			}
			edgenum++;
			if(u>vsize)
				vsize=u;
			if(v>vsize)
				vsize=v;

			edges.push_back(make_pair(u, v));
			if(!v_present[u]){
				_vsize++;
				v_present[u] = true;
			}
			if(!v_present[v]){
				_vsize++;
				v_present[v] = true;
			}
			if(u == 31147 && adj_graph.size() == 31148)
				cout<<"\t\tNormally something should be donne"<<endl;
		}
		i++;
	}
	int __vsize = u;
	vsize++;

	graph.resize(vsize+1);
	for(long long i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}
	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	
	sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{

			if(a.second<=b.second)
				return true;
			else
				return false;
		}

	});
	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}

	vector< pair<int, int> >().swap(edges);
	//edges.clear();
#ifndef Release	
	vector<int> inpos(vsize);
	for(int i=0; i<vsize; i++){
		inpos[i]=graph[i].instart;
	}
	inedge.resize(edgenum);
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outstart+graph[u].outdegree; j++){
			inedge[inpos[outedge[j]]]=u;
			inpos[outedge[j]]++;
		}
	}
#endif

	if(vsize!=adj_graph.size())
		cout <<"__vsize :"<<__vsize<<"_vsize :"<<_vsize<< "\tvsize: " << vsize <<"\t edgenum: "<<edgenum<<"\tInitial size: "<<adj_graph.size()<<"\tnb_edge: "<<nb_edge<<endl;
	/*cout << "vsize: " << vsize << endl;
	cout << "edgenum: " << edgenum <<"\tnb_edge: "<<nb_edge<<"\tedge.size():"<<edges.size()<< endl;*/
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;
}

void Graph_gorder::readGraph(const string& fullname) {
	FILE* fp;
	fp=fopen(fullname.c_str(), "r");
	if(fp==NULL){
		cout << "Fail to open " << fullname << endl;
		quit();
	}

	char line[40];
	int u, v;
	const char* str=NULL;
	
	vsize=0;
	edgenum=0;
	vector< pair<int, int> > edges;
	edges.reserve(100000000);

	while(feof(fp)!=true){
		if(fgets(line, 40, fp)){
			u=v=0;
			str=line;
			while(isdigit(*str))
				u=u*10+(*str++ - '0');
			str++;
			while(isdigit(*str))
				v=v*10+(*str++ - '0');

			if(u==v)
				continue;
			edgenum++;
			if(u>vsize)
				vsize=u;
			if(v>vsize)
				vsize=v;

			edges.push_back(make_pair(u, v));
		}
	}
	vsize++;

	fclose(fp);
	graph.resize(vsize+1);
	for(long long i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}
	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	
	sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{

			if(a.second<=b.second)
				return true;
			else
				return false;
		}

	});
	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}

	vector< pair<int, int> >().swap(edges);
#ifndef Release	
	vector<int> inpos(vsize);
	for(int i=0; i<vsize; i++){
		inpos[i]=graph[i].instart;
	}
	inedge.resize(edgenum);
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outstart+graph[u].outdegree; j++){
			inedge[inpos[outedge[j]]]=u;
			inpos[outedge[j]]++;
		}
	}
#endif

	cout << "vsize: " << vsize << endl;
	cout << "edgenum: " << edgenum << endl;
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;
}

void Graph_gorder::Transform(){
	vector<int> order;
	RCMOrder(order);
	if(order.size()!=vsize){
		cout << "order.size()!=vsize" << endl;
		quit();
	}
	if(graph.size()!=(vsize+1)){
		cout << "graph.size()!=(vsize+1)" << endl;
		quit();
	}

	vector<int>().swap(inedge);
	vector< pair<int, int> > edges;
	edges.reserve(edgenum);
	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart, limit=graph[i+1].outstart; j<limit; j++)
			edges.push_back(make_pair(order[i], order[outedge[j]]));
	}
	if(edges.size()!=edgenum){
		cout << "edges.size()!=edgenum" << endl;
		quit();
	}

	for(int i=0; i<vsize; i++){
		graph[i].outdegree=graph[i].indegree=0;
	}
	for(int i=0; i<edges.size(); i++){
		graph[edges[i].first].outdegree++;
		graph[edges[i].second].indegree++;
	}

	graph[0].outstart=0;
	graph[0].instart=0;
	for(int i=1; i<vsize; i++){
		graph[i].outstart=graph[i-1].outstart+graph[i-1].outdegree;
		graph[i].instart=graph[i-1].instart+graph[i-1].indegree;
	}
	graph[vsize].outstart=edgenum;
	graph[vsize].instart=edgenum;

	sort(edges.begin(), edges.end(), [](const pair<int, int>& a, const pair<int, int>& b)->bool{
		if(a.first<b.first)
			return true;
		else if(a.first>b.first)
			return false;
		else{

			if(a.second<=b.second)
				return true;
			else
				return false;
		}
	});

	outedge.resize(edgenum);
	for(long long i=0; i<edges.size(); i++){
		outedge[i]=edges[i].second;
	}
	vector< pair<int, int> >().swap(edges);
	vector<int> inpos(vsize);
	for(int i=0; i<vsize; i++){
		inpos[i]=graph[i].instart;
	}
	inedge.resize(edgenum);
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outstart+graph[u].outdegree; j++){
			inedge[inpos[outedge[j]]]=u;
			inpos[outedge[j]]++;
		}
	}
}


void Graph_gorder::writeGraph(ostream& out){
	for(int u=0; u<vsize; u++){
		for(int j=graph[u].outstart; j<graph[u].outdegree+graph[u].outstart; j++){
			int v=outedge[j];
			out << u << '\t' << v << endl;
		}
	}
}


vector<vector<int>> Graph_gorder::SetReOrderedGraph(vector<int>order){
	//vector<int>().swap(inedge);

	vector< vector<int> > ReOrderedGraph(vsize);
	int u, v, i;
  	//#pragma omp parallel for schedule(dynamic,1)
	for(int i=0; i<vsize; i++){
		u=order[i];
		ReOrderedGraph[u].reserve(graph[i+1].outstart-graph[i].outstart);
		for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
			v=order[outedge[j]];
			//cout<<"\t"<<v;
			ReOrderedGraph[u].push_back(v);
		}
		sort(ReOrderedGraph[u].begin(), ReOrderedGraph[u].end());
	}

/*	for(int u=0; u<ReOrderedGraph.size(); u++){
			cout << u << "(" <<ReOrderedGraph[u].size() << ") :";
		for(int j=0; j<ReOrderedGraph[u].size(); j++){
			cout << ReOrderedGraph[u][j] <<"\t";
		}
		cout<<endl;
	}*/
	//#pragma omp barrier
	return ReOrderedGraph;
}

void Graph_gorder::PrintReOrderedGraph(const vector<int>& order){
	ofstream out((name+"_Gorder.txt").c_str());

	vector<int>().swap(inedge);

	vector< vector<int> > ReOrderedGraph(vsize);
	int u, v;
	for(int i=0; i<vsize; i++){
		u=order[i];
		ReOrderedGraph[u].reserve(graph[i+1].outstart-graph[i].outstart);
		for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
			v=order[outedge[j]];
			ReOrderedGraph[u].push_back(v);
		}
		sort(ReOrderedGraph[u].begin(), ReOrderedGraph[u].end());
	}
/*
	for(int u=0; u<vsize; u++){
		sort(ReOrderedGraph[u].begin(), ReOrderedGraph[u].end());
	}
*/
	for(int u=0; u<vsize; u++){
		for(int j=0; j<ReOrderedGraph[u].size(); j++){
			out << u << '\t' << ReOrderedGraph[u][j] << endl;
		}
	}
	out.close();
}


void Graph_gorder::GraphAnalysis(){
	vector<int> tmp(vsize);
	for(int i=0; i<vsize; i++){
		tmp[i]=graph[i].outdegree;
	}
	sort(tmp.begin(), tmp.end());
	
	cout << "outdegree:" << endl;
	vector<int>::iterator tmpit1, tmpit2;
	tmpit1=tmp.begin();
	for(int i=1; i<vsize; i*=10){
		tmpit2=tmpit1;
		tmpit1=upper_bound(tmp.begin(), tmp.end(), i);
		cout << i << ": " << tmpit1-tmpit2 << endl;
	}


	for(int i=0; i<vsize; i++){
		tmp[i]=graph[i].indegree;
	}
	sort(tmp.begin(), tmp.end());
	
	cout << "indegree:" << endl;
	tmpit1=tmp.begin();
	for(int i=1; i<vsize; i*=10){
		tmpit2=tmpit1;
		tmpit1=upper_bound(tmp.begin(), tmp.end(), i);
		cout << i << ": " << tmpit1-tmpit2 << endl;
	}
}


void Graph_gorder::RemoveDuplicate(const string& fullname) {
	FILE* fp;
	fp=fopen(fullname.c_str(), "r");
	if(fp==NULL){
		cout << "Fail to open " << fullname << endl;
		quit();
	}

	char line[40];
	int u, v;
	const char* str=NULL;
	
	set< pair<int, int> > edges;

	while(feof(fp)!=true){
		if(fgets(line, 40, fp)){
			u=v=0;
			str=line;
			while(isdigit(*str))
				u=u*10+(*str++ - '0');
			str++;
			while(isdigit(*str))
				v=v*10+(*str++ - '0');

			if(u==v)
				continue;

			edges.insert(make_pair(u, v));
		}
	}
	
	fclose(fp);

	cout << "after remove, the size is " << edges.size() << endl;
	
	ofstream fout;
	fout.open("NoDuplicate.txt");
	for(set< pair<int, int> >::iterator it=edges.begin(); it!=edges.end(); it++){
		fout << it->first << '\t' << it->second << endl;
	}
	fout.close();
}


vector<string>& Graph_gorder::split(const string &s, char delim, vector<string> &elems) {
	int begin, end;

	begin=0;
	end=s.find(delim);
	while(end!=string::npos){
		elems.push_back(s.substr(begin, end-begin));
		begin=end+1;
		end=s.find(delim, begin);
	}
	if(begin!=s.size()){
		elems.push_back(s.substr(begin));
	}

	return elems;
}

vector<string> Graph_gorder::split(const string &s, char delim) {
	vector<string> elems;
	return split(s, delim, elems);
}

void Graph_gorder::strTrimRight(string& str) {
	string whitespaces(" \t\r\n");
	int index = str.find_last_not_of(whitespaces);
	if (index != string::npos) 
		str.erase(index+1);
	else
		str.clear();
}

void Graph_gorder::GapCount(){
	int* gap = new int[vsize];
	memset(gap, 0, sizeof(int)*vsize);

	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart+1; j<graph[i].outdegree+graph[i].outstart; j++){
				gap[outedge[j]-outedge[j-1]]++;
		}
		gap[outedge[graph[i].outstart]]++;
	}

	double entropy=0;
	for(int i=0; i<vsize; i++){
		if(gap[i]==0)
			continue;
		else{
			entropy+=(double)gap[i]/edgenum*log2((double)gap[i]/edgenum);
		}
	}
	cout << "shannon: " << entropy << endl;
	delete[] gap;
}

double Graph_gorder::GapCost(vector<int>& order){
	double gaplog=0;
	double gaplog2=0;
	vector<int> edgelist;
	edgelist.reserve(100000);
	for(int i=0; i<vsize; i++){
		for(int j=graph[i].outstart+1; j<graph[i].outdegree+graph[i].outstart; j++){
			if(outedge[j]-outedge[j-1])
				gaplog+=log(double(outedge[j]-outedge[j-1]))/log(double(2));
		}
		edgelist.clear();
		for(int j=graph[i].outstart; j<graph[i].outstart+graph[i].outdegree; j++){
			edgelist.push_back(order[outedge[j]]);
		}
		sort(edgelist.begin(), edgelist.end());
		for(int j=1; j<edgelist.size(); j++){
			if(edgelist[j]-edgelist[j-1])
				gaplog2+=log(double(edgelist[j]-edgelist[j-1]))/log(double(2));
		}
	}
	cout << "original average gap cost: " << gaplog/edgenum << endl;
	cout << "new average gap cost: " << gaplog2/edgenum << endl;

	return gaplog2/edgenum;
}


void Graph_gorder::GorderGreedy(vector<int>& retorder, int window){
	UnitHeap unitheap(vsize);
	vector<bool> popvexist(vsize, false);
	vector<int> order;
	int count=0;
	vector<int> zero;
	zero.reserve(10000);
	order.reserve(vsize);
	const int hugevertex=sqrt((double)vsize);

	for(int i=0; i<vsize; i++){
		unitheap.LinkedList[i].key=graph[i].indegree;
		unitheap.update[i]=-graph[i].indegree;
	}
	unitheap.ReConstruct();

	int tmpindex, tmpweight;

	tmpweight=-1;
	for(int i=0; i<vsize; i++){
		if(graph[i].indegree>tmpweight){
			tmpweight=graph[i].indegree;
			tmpindex=i;
		}else if(graph[i].indegree+graph[i].outdegree==0){
			unitheap.update[i]=INT_MAX/2;
			zero.push_back(i);
			unitheap.DeleteElement(i);
		}
	}

	order.push_back(tmpindex);
	unitheap.update[tmpindex]=INT_MAX/2;
	unitheap.DeleteElement(tmpindex);
	for(int i=graph[tmpindex].instart, limit1=graph[tmpindex+1].instart; i<limit1; i++){
		int u=inedge[i];
		if(graph[u].outdegree<=hugevertex){
			if(unitheap.update[u]==0){
				unitheap.IncrementKey(u);
			} else {
#ifndef Release
				if(unitheap.update[u]==INT_MAX)
					unitheap.update[u]=INT_MAX/2;
#endif
				unitheap.update[u]++;
			}
			
			if(graph[u].outdegree>1)
			for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
				int w=outedge[j];
				if(unitheap.update[w]==0){
					unitheap.IncrementKey(w);
				} else {
#ifndef Release
					if(unitheap.update[w]==INT_MAX)
						unitheap.update[w]=INT_MAX/2;
#endif
					unitheap.update[w]++;
				}
				
			}
		}
	}
	if(graph[tmpindex].outdegree<=hugevertex){
		for(int i=graph[tmpindex].outstart, limit1=graph[tmpindex+1].outstart; i<limit1; i++){
			int w=outedge[i];
			if(unitheap.update[w]==0){
				unitheap.IncrementKey(w);
			}else{
#ifndef Release
				if(unitheap.update[w]==INT_MAX)
					unitheap.update[w]=INT_MAX/2;
#endif
				unitheap.update[w]++;
			}
			
		}
	}

#ifndef Release
	clock_t time1, time2, time3, time4;
	clock_t sum1=0, sum2=0, sum3=0;
#endif

	while(count<vsize-1-zero.size()){
#ifndef Release
		if(count%1000000==0){
			/*cout << count << endl;
			cout << "sum1: " << sum1 << endl;
			cout << "sum2: " << sum2 << endl;
			cout << "sum3: " << sum3 << endl;
			cout << endl;*/
			sum1=sum2=sum3=0;
		}
				
		time1=clock();
#endif

		int v=unitheap.ExtractMax();
		count++;
#ifndef Release
		time2=clock();
#endif
		order.push_back(v);
		unitheap.update[v]=INT_MAX/2;

		int popv;
		if(count-window>=0)
			popv=order[count-window];
		else
			popv=-1;

		if(popv>=0){
			if(graph[popv].outdegree<=hugevertex){
				for(int i=graph[popv].outstart, limit1=graph[popv+1].outstart; i<limit1; i++){
					int w=outedge[i];
					unitheap.update[w]--;
#ifndef Release
					if(unitheap.update[w]==0)
						unitheap.update[w]=INT_MAX/2;
#endif
				}
			}

			for(int i=graph[popv].instart, limit1=graph[popv+1].instart; i<limit1; i++){
				int u=inedge[i];
				if(graph[u].outdegree<=hugevertex){
					unitheap.update[u]--;
#ifndef Release
					if(unitheap.update[u]==0)
						unitheap.update[u]=INT_MAX/2;
#endif
					if(graph[u].outdegree>1)
					if(binary_search(&outedge[graph[u].outstart], &outedge[graph[u+1].outstart], v)==false){
						for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
							int w=outedge[j];
							unitheap.update[w]--;
#ifndef Release
							if(unitheap.update[w]==0)
								unitheap.update[w]=INT_MAX/2;
#endif
						}
					} else {
						popvexist[u]=true;
					}
				}
			}
		}

#ifndef Release
		time3=clock();
#endif
		if(graph[v].outdegree<=hugevertex){
			for(int i=graph[v].outstart, limit1=graph[v+1].outstart; i<limit1; i++){
				int w=outedge[i];
				if(__builtin_expect(unitheap.update[w]==0, 0)){
					unitheap.IncrementKey(w);
				} else {
#ifndef Release
					if(unitheap.update[w]==INT_MAX)
						unitheap.update[w]=INT_MAX/2;
#endif
					unitheap.update[w]++;
				}
				
			}
		}

		for(int i=graph[v].instart, limit1=graph[v+1].instart; i<limit1; i++){
			int u=inedge[i];
			if(graph[u].outdegree<=hugevertex){
				if(__builtin_expect(unitheap.update[u]==0, 0)){
					unitheap.IncrementKey(u);
				} else {
#ifndef Release
					if(unitheap.update[u]==INT_MAX)
						unitheap.update[u]=INT_MAX/2;
#endif
					unitheap.update[u]++;
				}
				
				if(popvexist[u]==false){
					if(graph[u].outdegree>1)
					for(int j=graph[u].outstart, limit2=graph[u+1].outstart; j<limit2; j++){
						int w=outedge[j];
						if(__builtin_expect(unitheap.update[w]==0, 0)){
							unitheap.IncrementKey(w);
						}else{
#ifndef Release
							if(unitheap.update[w]==INT_MAX)
								unitheap.update[w]=INT_MAX/2;
#endif
							unitheap.update[w]++;
						}
					}
				} else {
					popvexist[u]=false;
				}
			}
		}


#ifndef Release
	time4=clock();
	sum1+=time2-time1;
	sum2+=time3-time2;
	sum3+=time4-time3;
#endif
	}
	order.insert(order.end()-1, zero.begin(), zero.end());


#ifndef Release
	vector<int> tmporder=order;
	sort(tmporder.begin(), tmporder.end());
	for(int i=0; i<tmporder.size()-1; i++){
		if(tmporder[i]==tmporder[i+1]){
			cout << "same elements: " << tmporder[i] << endl;
			system("pause");
		}
	}
	for(int i=0; i<tmporder.size(); i++){
		if(tmporder[i]!=i){
			cout << tmporder[i] << '\t' << i << endl;
			system("pause");
		}
	}
	tmporder=vector<int>();
#endif

	retorder.clear();
	retorder.resize(vsize);
	for(int i=0; i<vsize; i++){
		retorder[order[i]]=i;
	}
}


void Graph_gorder::RCMOrder(vector<int>& retorder){
	queue<int> que;
	bool* BFSflag=new bool[vsize];
	bool* QueFlag=new bool[vsize];
	memset(BFSflag, 0, sizeof(bool)*vsize);
	memset(QueFlag, 0, sizeof(bool)*vsize);

	vector<int> tmp;
	vector<int> degreevertex(vsize);
	for(int i=0; i<vsize; i++){
		degreevertex[i]=i;
	}

	sort(degreevertex.begin(), degreevertex.end(), [&](const int& a, const int& b)->bool{
		if(graph[a].outdegree+graph[a].indegree<graph[b].outdegree+graph[b].indegree)
			return true;
		else
			return false;
	});

        int now;
	vector<int> order;

	for(int k=0; k<vsize; k++){
		int i=degreevertex[k];
		if(BFSflag[i]==false){
			que.push(i);
//			QueFlag[i]=true;
			BFSflag[i]=true;
			order.push_back(i);

			while(que.empty()==false){
				now=que.front();
				que.pop();

//				BFSflag[now]=true;
				tmp.clear();
				if(now >= graph.size())
				{
					cout << "I wander why it >= graph.size, I should stop"<<endl;
					//quit();
				}
				for(int it=graph[now].outstart, limit=graph[now+1].outstart; it<limit; it++){
					tmp.push_back(outedge[it]);
					//cout << "\t it = "<< it << "\tinserting "<<outedge[it]<< "\tlimit = " << limit<< endl;
					if(outedge[it] >= graph.size())
					{
						cout << "graph.size()= "<< graph.size()<<" outedge.size() =" << outedge.size()<<"\n  now="<<now<<" graph["<<now<<"].outstart=" <<graph[now].outstart<< endl;
						cout << "now+1="<<(now+1)<<" graph["<<(now+1)<<"].outstart=" <<graph[now+1].outstart<< endl;
						cout << "I wander why it >= graph.size, I should stop"<<endl;
						//quit();
					}
				}
				sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b)->bool{
					if(graph[a].outdegree+graph[a].indegree<graph[b].outdegree+graph[b].indegree)
						return true;
					else
						return false;
				});
				if(tmp.size()!=graph[now].outdegree)
				{
					cout << "tmp.size()="<<tmp.size()<<"!=graph["<<now<<"].outdegree=" <<graph[now].outdegree<< endl;

				}

				for(int i=0; i<tmp.size(); i++){
//					if((BFSflag[tmp[i]]==false)&&(QueFlag[tmp[i]]==false)){
					if(BFSflag[tmp[i]]==false){
						que.push(tmp[i]);
                        			BFSflag[tmp[i]]=true;
						order.push_back(tmp[i]);
                        		}
					//cout << "\t Last for in  while "<<i<< endl;
				} 
				//cout << "dans while "<<que.size()<< endl;
        		}
		}
	}

        delete[] BFSflag;
        delete[] QueFlag;

	if(order.size()!=vsize){
		cout << "order.size()!=vsize" << endl;
		quit();
	}

	retorder.resize(vsize);
	for(int i=0; i<order.size(); i++){
		retorder[order[i]]=order.size()-1-i;
	}
}


unsigned long long Graph_gorder::LocalityScore(const int w){
	unsigned long long sum=0;
	for(int i=0; i<vsize; i++){
		for(int j=i-1; j>=i-w && j>=0; j--){
			sum+=IntersectionSize(inedge.data()+graph[i].instart, inedge.data()+graph[j].instart, graph[i].indegree, graph[j].indegree, -1);
			if(binary_search(inedge.data()+graph[i].instart, inedge.data()+graph[i].instart+graph[i].indegree, j))
				sum++;
			if(binary_search(inedge.data()+graph[j].instart, inedge.data()+graph[j].instart+graph[j].indegree, i))
				sum++;
		}
	}
	return sum;
}

