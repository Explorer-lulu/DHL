#include "head.h"
#include<stdio.h>
//#include<metis.h>
#include<vector>
#include<stdlib.h>
#include<memory.h>
#include<unordered_map>
#include<map>
#include<set>
#include<deque>
#include<stack>
#include<algorithm>
#include<sys/time.h>

// HGP-tree build

#define PARTITION_PART_NF 4
#define LEAF_TAU 16
#define FILE_ONTREE_MIND	  "fill your file name here"

typedef int vertex_id;
int bufrSize;  //Setting
int siInt = sizeof(int);
int intPrBlo;  // e.g., bufrSize/siInt;
int nodtid;

typedef struct{
	double x,y;
	vector<int> adjnodes;
	vector<int> adjweight;
	bool isborder;
	vector<int> HGPtreepath; 
}HGPNode;

typedef struct{
	vector<int> borders;
	vector<int> children;
	bool isleaf;
	vector<int> leafnodes;
	int father;
	vector<int> union_borders; 
	vector<int> mind; 
	vector<int> nonleafinvlist;
	vector<int> leafinvlist;
	vector<int> up_pos;
	vector<int> current_pos;
}HGPTreeNode;

class WriteBuffer {
public:
	FILE *outfile;
	int *buf;
	int ptr;

public:
	WriteBuffer() {
	};
	WriteBuffer(FILE *outfile_t) {
		buf = (int*)malloc(bufrSize);
		open(outfile_t);
	};
	~WriteBuffer() {
		flush();
		free(buf);
	};

	void write(int *src, int cnt) {
		unsigned total=cnt;
		while(cnt>0){
		        if (cnt < intPrBlo-ptr) {
		                memcpy(buf+ptr,src+(total-cnt),cnt*siInt);
		                ptr += cnt;
		                cnt = 0;
		        } else {
		                memcpy(buf+ptr,src+(total-cnt),(intPrBlo-ptr)*siInt);
		                fwrite(buf,siInt,intPrBlo,outfile);
		                cnt -= intPrBlo-ptr;
		                ptr = 0;
		        }
		}
	};

	void open(FILE *outfile_t) {
		outfile = outfile_t;
		ptr = 0;
		rewind(outfile);
	}

	void flush() {
		if(ptr>0) {
		        fwrite(buf,siInt,ptr,outfile);
		        ptr = 0;
		}
	};
};

int noe; 
vector<HGPNode> HGPNodes;
vector<HGPTreeNode> HGPTree;

typedef struct{
	int tnid; 
	set<int> nset; 
}Status;

set<int> childset[PARTITION_PART_NF];

void Graph::HGPbuild(){

	HGPTreeNode root;
	root.isleaf = false;
	root.father = -1;
	HGPTree.push_back(root);

	stack<Status> buildstack;
	Status rootstatus;
	rootstatus.tnid = 0;
	rootstatus.nset.clear();
	for ( int i = 0; i < HGPNodes.size(); i++ ){
		rootstatus.nset.insert(i);
	}
	buildstack.push( rootstatus );

	// start to build
	unordered_map<int,int> presult;
	set<int> childset[PARTITION_PART_NF];


	while( buildstack.size() > 0 ){
		Status current = buildstack.top();
		buildstack.pop();

		for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
			HGPNodes[*it].HGPtreepath.push_back( current.tnid );
		}
		if ( current.nset.size() <= LEAF_TAU ){
			HGPTree[current.tnid].isleaf = true;
			HGPTree[current.tnid].leafnodes.clear();
			for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
				HGPTree[current.tnid].leafnodes.push_back( *it );
			}
			continue;
		}
		for ( int i = 0; i < PARTITION_PART_NF; i++ ){
			childset[i].clear();
		}
		int slot;
		for ( set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
			slot = presult[*it];
			childset[slot].insert(*it);
		}

		int childpos;
		for ( int i = 0; i < PARTITION_PART_NF; i++ ){
			HGPTreeNode tnode;
			tnode.isleaf = false;
			tnode.father = current.tnid;
			

            HGPTree.push_back(tnode);
			childpos = HGPTree.size() - 1;
			HGPTree[current.tnid].children.push_back( childpos );


			HGPTree[childpos].borders.clear();
			for ( set<int>::iterator it = childset[i].begin(); it != childset[i].end(); it++ ){

				bool isborder = false;
				for ( int j = 0; j < HGPNodes[*it].adjnodes.size(); j++ ){
					if ( childset[i].find( HGPNodes[*it].adjnodes[j] ) == childset[i].end() ){
						isborder = true;
						break;
					}
				}
				if ( isborder ){
					HGPTree[childpos].borders.push_back(*it);
					HGPNodes[*it].isborder = true;
				}
			}

			Status ongoingstatus;
			ongoingstatus.tnid = childpos;
			ongoingstatus.nset = childset[i];
			buildstack.push(ongoingstatus);

		}

	}
}

void Graph::hierarchy_shortestd_for_node_load(){
	FILE* filein = fopen( FILE_ONTREE_MIND, "rb" );
	int* buf;
	int count, pos = 0;
	while( fread( &count, sizeof(int), 1, filein ) == 1){

		buf = new int[count];
		if(fread( buf, sizeof(int), count, filein )==1);
		HGPTree[pos].union_borders.clear();
		for ( int i = 0; i < count; i++ ){
			HGPTree[pos].union_borders.push_back(buf[i]);
		}
		delete[] buf;
		if(fread( &count, sizeof(int), 1, filein )==1);
		buf = new int[count];
		if(fread( buf, sizeof(int), count, filein )==1);
		HGPTree[pos].mind.clear();
		for ( int i = 0; i < count; i++ ){
			HGPTree[pos].mind.push_back(buf[i]);
		}
		pos++;
		delete[] buf;
	}
	fclose(filein);
}

// DHL

struct weightedge{
	vertex_id nodtid;
	int weight;
};

void labelInit(int &u, int &u_outdeg, int &u_indeg, weightedge * neighbor, WriteBuffer & write_label, int *buf_out, int *buf_in)
{
	int ptr_out = 0, ptr_in = 0;
	buf_out[ptr_out++] = u;
	buf_out[ptr_out++] = u_outdeg+1;
	buf_out[ptr_out++] = u_indeg+1;
	bool finish = false;
	for(int i = 0; i < u_outdeg; i++)
	{
		if(!finish)
		{
			if(neighbor[i].nodtid > u)
			{
				finish = true;
				buf_out[ptr_out++] = u;
				buf_out[ptr_out++] = 0;
			}
		}
		buf_out[ptr_out++] = neighbor[i].nodtid;
		buf_out[ptr_out++] = neighbor[i].weight;
	}
	if(!finish)
	{
		finish = true;
		buf_out[ptr_out++] = u;
		buf_out[ptr_out++] = 0;
	}
	finish = false;
	for(int i=u_outdeg; i<u_indeg+u_outdeg; i++){
		if(!finish)
		{
			if(neighbor[i].nodtid > u)
			{
				finish = true;
				buf_in[ptr_in++] = u;
				buf_in[ptr_in++] = 0;
			}
		}
		buf_in[ptr_in++] = neighbor[i].nodtid;
		buf_in[ptr_in++] =neighbor[i].weight;
	}
	if(!finish)
	{
		finish = true;
		buf_in[ptr_in++] = u;
		buf_in[ptr_in++] = 0;
	}
	write_label.write(buf_out, 3+2*(u_outdeg+1));
	write_label.write(buf_in, 2*(u_indeg+1));
	write_label.flush();
}

vector<int> _DD,_DD2;
struct DegComp{
	int x;
	DegComp(int _x){
		x=_x;
	}
	bool operator< (const DegComp d) const{
		if(_DD[x]!=_DD[d.x])
			return _DD[x]<_DD[d.x];
		if(_DD2[x]!=_DD2[d.x])
			return _DD2[x]<_DD2[d.x];
		return x<d.x;
	}
};


vector<int> NodeOrders;
struct OrderComp{
	int x;
	int y;
	OrderComp(int _x, int _y){
		x=_x; y=_y;
	}
	bool operator< (const OrderComp& d) const{
		if(x==d.x && y==d.y){
			return false;
		}else{
			if(x!=d.x)
				return NodeOrders[x]<NodeOrders[d.x];
			if(y!=d.y)
				return NodeOrders[y]<NodeOrders[d.y];
		}
	}
};

void Graph::DHLGenerate(){
	int Twidth=0;

	map<int, vector<int>> mi;
	SCconNodesMT.assign(nodenum, mi);

	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbor.size();i++){
		for(int j=0;j<Neighbor[i].size();j++)
			E[i].insert(make_pair(Neighbor[i][j].first,make_pair(0,1)));
	}

	_DD.assign(nodenum,0);_DD2.assign(nodenum,0);
	DD.assign(nodenum,0); DD2.assign(nodenum,0);

	set<DegComp> Deg;
	int degree;
	for(int i=0;i<nodenum;i++){
		degree=Neighbor[i].size();
		if(degree!=0){
			_DD[i]=degree;
			_DD2[i]=degree;
			DD[i]=degree;
			DD2[i]=degree;
			Deg.insert(DegComp(i));
		}
	}

	vector<bool> exist; exist.assign(nodenum,true);
	vector<bool> change; change.assign(nodenum,false);

	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);

	cout<<"Begin to contract vertex"<<endl;
	int count=0;

	while(!Deg.empty()){
		if(count%10000==0)
			cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
		count+=1;
		int x=(*Deg.begin()).x;

		while(true){
			if(change[x]){
				Deg.erase(DegComp(x));
				_DD[x]=DD[x];
				_DD2[x]=DD2[x];
				Deg.insert(DegComp(x));
				change[x]=false;
				x=(*Deg.begin()).x;
			}else
				break;
		}

		vNodeOrder.push_back(x);
		Deg.erase(Deg.begin());
		exist[x]=false;

		vector<pair<int,pair<int,int>>> Neigh; 

		for(auto it=E[x].begin();it!=E[x].end();it++){
			if(exist[(*it).first]){
				Neigh.push_back(*it);
			}
		}

		if(Neigh.size()>Twidth)
			Twidth=Neigh.size();

		NeighborCon[x].assign(Neigh.begin(),Neigh.end());

		for(int i=0;i<Neigh.size();i++){
			int y=Neigh[i].first;
			deleteE(x,y);
			change[y]=true;
		}

		int stepf=Neigh.size()/threadnum;
		boost::thread_group threadf;
		for(int i=0;i<threadnum;i++){
			pair<int,int> p;
			p.first=i*stepf;
			if(i==threadnum-1)
				p.second=Neigh.size();
			else
				p.second=(i+1)*stepf;
			threadf.add_thread(new boost::thread(&Graph::NeighborComOrderGenerate, this, boost::ref(Neigh), p, x));
		}
		threadf.join_all();
	}

	NodeOrder.assign(nodenum,-1);
	for(int k=0;k<vNodeOrder.size();k++){
		NodeOrder[vNodeOrder[k]]=k;
	}
	cout<<"Finish Contract"<<" , treewidth "<<Twidth<<endl;
}


int Graph::StartDHL(string indexfile, string orderfile){
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	fstream file;
	file.open(indexfile);
	if(!file){
		t1=std::chrono::high_resolution_clock::now();
		DHLconsorder(orderfile);
		t2=std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
		runT= time_span.count();
		cout<<"DHL order MT "<<runT<<endl;
		writeShortCutDouble(indexfile);
		ReadShortCutDouble(indexfile);
	}
	else
	{
		ReadShortCutDouble(indexfile);
	}

	return 0;
}

void Graph::insertEdge(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
	}
	else{
		if(E[u][v].first>w)
			E[u][v]=make_pair(w,1);
		else if(E[u][v].first==w)
			E[u][v].second+=1;
	}
}

void Graph::DHLconsorder(string orderfile){
	ifstream IF(orderfile);
	if(!IF){
		cout<<"Cannot open Map "<<orderfile<<endl;
	}
	NodeOrder.assign(nodenum, -1);
	vNodeOrder.assign(nodenum, -1);
	int num, nodeID, nodeorder;
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>nodeID>>nodeorder;
		NodeOrder[nodeID]=nodeorder;
		if(nodeorder!=-1){
			vNodeOrder[nodeorder]=nodeID;
		}
	}

	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);

	map<int, vector<int>> mi;
	SCconNodesMT.assign(nodenum, mi);

	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbor.size();i++){
		for(int j=0;j<Neighbor[i].size();j++)
			E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
	}

	vector<bool> exist; exist.assign(nodenum,true);

	for(int nodeorder=0;nodeorder<nodenum;nodeorder++){
		int x=vNodeOrder[nodeorder];
		if(x!=-1){
			exist[x]=false;

			vector<pair<int,pair<int,int>>> Neigh; 

			for(auto it=E[x].begin();it!=E[x].end();it++){
				if(exist[(*it).first]){
					Neigh.push_back(*it);
				}
			}
			NeighborCon[x].assign(Neigh.begin(),Neigh.end());

			for(int i=0;i<Neigh.size();i++){
				int y=Neigh[i].first;
				deleteEorder(x,y);
			}

			if(Neigh.size()<=100){
				for(int i=0;i<Neigh.size();i++){
					for(int j=i+1;j<Neigh.size();j++){
						insertEorder(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
						if(Neigh[i].first<Neigh[j].first)
							SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);
						else if(Neigh[j].first<Neigh[i].first)
							SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);
					}
				}
			}else{
				int step=Neigh.size()/threadnum;
				boost::thread_group thread;
				for(int i=0;i<threadnum;i++){
					pair<int,int> p;
					p.first=i*step;
					if(i==threadnum-1)
						p.second=Neigh.size();
					else
						p.second=(i+1)*step;
					thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
				}
				thread.join_all();
			}

		}
	}
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
	sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		w1=Neighvec[k].second.first;
		for(int h=0;h<Neighvec.size();h++){
			ID2=Neighvec[h].first;
			w2=Neighvec[h].second.first;
			insertEdge(ID1, ID2, w1+w2);
			if(ID1<ID2)
				SCconNodesMT[ID1][ID2].push_back(x);
		}
	}
	sm->notify();
}


int Graph::DHLcontractionlevel(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult){
	benchmark::heap<2,int,int> Heap(nodenum);
	vector<int> vDistance(nodenum, INF);
	int topNodeID, topDistance, neighborNodeID, neighborLength;
	int w;
	vDistance[ID1]=0;
	Heap.update(ID1,0);
	map<int,int> mWDistance;
	vector<pair<int,int>> vpWDistance;
	map<int,int> mDistance;
	int maxWDistance=-1;

	for(auto ivp=vW.begin();ivp!=vW.end();ivp++){
		w=(*ivp).first;
		int d=(*ivp).second+dUV;
		mWDistance[w]=d;
		if(d>maxWDistance)
			maxWDistance=d;
		mDistance[w]=INF;
	}

	if(mDistance.empty()){
		return 0;
	}

	for(auto imWDistance=mWDistance.begin();imWDistance!=mWDistance.end();imWDistance++)
		vpWDistance.push_back(make_pair((*imWDistance).first, (*imWDistance).second));

	int dThreshold=maxWDistance;

	while(!Heap.empty()){
		Heap.extract_min(topNodeID, topDistance);
		if(vbVisited[topNodeID])
			continue;
		if(topDistance>dThreshold)
			break;
		for(auto ivp=vvNode[topNodeID].begin();ivp!=vvNode[topNodeID].end();ivp++){
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second;
			if(vbVisited[neighborNodeID] || neighborNodeID==ID2)
				continue;
			int d = vDistance[topNodeID] + neighborLength;
			if(vDistance[neighborNodeID] ==INF && neighborNodeID != ID1)
			{
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}else if(vDistance[neighborNodeID]>d){
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
			else if(vDistance[neighborNodeID] == d)
			{
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
		}
	}

	for(auto imDistance = mDistance.begin(); imDistance != mDistance.end(); imDistance++)
	{
		if((*imDistance).second > mWDistance[(*imDistance).first])
		{
			int w = (*imDistance).first;
			int distance = mWDistance[(*imDistance).first];

			vvpResult[k].push_back(make_pair(w, distance));
		}
	}

	return 0;
}


int Graph::DHLsyncWP(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult){
	benchmark::heap<2,int,int> Heap(nodenum);
	vector<int> vDistance(nodenum, INF);
	vector<int> vPre(nodenum,-1);
	int topNodeID, topDistance, neighborNodeID, neighborLength;
	int w;
	vDistance[ID1]=0;
	Heap.update(ID1,0);
	map<int,int> mWDistance;
	vector<pair<int,int>> vpWDistance;
	map<int,int> mDistance;
	int maxWDistance=-1;

	for(auto ivp=vW.begin();ivp!=vW.end();ivp++){
		w=(*ivp).first;
		int d=(*ivp).second+dUV;
		mWDistance[w]=d;
		if(d>maxWDistance)
			maxWDistance=d;
		mDistance[w]=INF;
	}

	if(mDistance.empty()){
		return 0;
	}

	for(auto imWDistance=mWDistance.begin();imWDistance!=mWDistance.end();imWDistance++)
		vpWDistance.push_back(make_pair((*imWDistance).first, (*imWDistance).second));

	int dThreshold=maxWDistance;

	while(!Heap.empty()){
		Heap.extract_min(topNodeID, topDistance);
		if(vbVisited[topNodeID])
			continue;
		if(topDistance>dThreshold)
			break;
		for(auto ivp=vvNode[topNodeID].begin();ivp!=vvNode[topNodeID].end();ivp++){
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second;
			if(vbVisited[neighborNodeID] || neighborNodeID==ID2)
				continue;
			int d = vDistance[topNodeID] + neighborLength;
			if(vDistance[neighborNodeID] ==INF && neighborNodeID != ID1)
			{
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				vPre[neighborNodeID]=topNodeID;
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}else if(vDistance[neighborNodeID]>d){
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				vPre[neighborNodeID]=topNodeID;
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
			else if(vDistance[neighborNodeID] == d)
			{
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
		}
	}


	for(auto imDistance = mDistance.begin(); imDistance != mDistance.end(); imDistance++)
	{
		int w = (*imDistance).first;
		if((*imDistance).second > mWDistance[(*imDistance).first])
		{

			int distance = mWDistance[(*imDistance).first];

			vvpResult[k].push_back(make_pair(w, distance));
		}
		else{

			vector<int> path;
			path.push_back(w);
			int pre=vPre[w];

			while(pre!=ID1){
				path.insert(path.begin(), pre);
				pre=vPre[pre];
			}
			path.insert(path.begin(), pre);


			for(int i=0;i<path.size()-1;i++){
				int currentID = path[i];
				while(OutEdgesM[currentID].find(path[i+1])!=OutEdgesM[currentID].end()){
					int MnodeID=OutEdgesM[currentID][path[i+1]];
					path.insert(path.begin()+i+1, MnodeID);
				}
			}
			int a,b;
			for(int k=0;k<path.size()-1;k++){
				if(path[k]<path[k+1]){
					a=path[k];
					b=path[k+1];
				}else{
					a=path[k+1];
					b=path[k];
				}
				int edgeID=EdgeRe[make_pair(a,b)];
				vSm[edgeID]->wait();
				tri TRI;
				TRI.u=ID1;
				TRI.v=ID2;
				TRI.w=w;
				EdgeOnPath[edgeID].push_back(TRI);
				vSm[edgeID]->notify();
			}

			PathInfor[ID1].insert(make_pair(make_pair(ID2,w),(*imDistance).second));

		}
	}

	return 0;
}


void Graph::insertEorder(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
	}
	else{
		if(E[u][v].first>w)
			E[u][v]=make_pair(w,1);
		else if(E[u][v].first==w)
			E[u][v].second+=1;
	}

	if(E[v].find(u)==E[v].end()){
		E[v].insert(make_pair(u,make_pair(w,1)));
	}
	else{
		if(E[v][u].first>w)
			E[v][u]=make_pair(w,1);
		else if(E[v][u].first==w)
			E[v][u].second+=1;
	}
}

void Graph::deleteEorder(int u,int v){
	if(E[u].find(v)!=E[u].end()){
		E[u].erase(E[u].find(v));
	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
	}
}



void Graph::deleteE(int u,int v){
	if(E[u].find(v)!=E[u].end()){
		E[u].erase(E[u].find(v));
		DD[u]--;
	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
		DD[v]--;
	}
}

int Graph::NewSCweight(int s, int t){
	int wt=INF;
	for(int i=0;i<Neighbor[s].size();i++){
		if(Neighbor[s][i].first==t){
			wt=Neighbor[s][i].second;
			break;
		}
	}
	int ssw,wtt,wid;
	vector<int> Wnodes; 
	if(s<t)
		Wnodes=SupportNodes[s][t];
	else
		Wnodes=SupportNodes[t][s];
	for(int i=0;i<Wnodes.size();i++){
		wid=Wnodes[i];
		for(int j=0;j<AdjaShort[wid].size();j++){
			if(AdjaShort[wid][j].first==s){
				ssw=AdjaShort[wid][j].second;
			}
			if(AdjaShort[wid][j].first==t){
				wtt=AdjaShort[wid][j].second;
			}
		}
		if(ssw+wtt<wt){
			wt=ssw+wtt;
		}
	}
	return wt;
}

int Graph::DHLcontractInc(int ID1, int ID2, int dUV, vector<pair<int, int> >& vW, vector<pair<int,int>>& addedShortcut){
	benchmark::heap<2,int,int> Heap(nodenum);
	vector<int> vDistance(nodenum, INF);
	int topNodeID, topDistance, neighborNodeID, neighborLength;
	int w;
	vDistance[ID1]=0;
	Heap.update(ID1,0);
	int count=0;
	map<int,int> mWDistance;
	vector<pair<int,int>> vpWDistance;
	map<int,int> mDistance;
	int maxWDistance=-1;

	for(auto ivp=vW.begin();ivp!=vW.end();ivp++){
		w=(*ivp).first;
		if(w!=ID1){
			int d=(*ivp).second+dUV;
			mWDistance[w]=d;
			if(d>maxWDistance)
				maxWDistance=d;
			mDistance[w]=INF;
		}
	}
	if(mDistance.empty()){
		return 0;
	}
	for(auto imWDistance=mWDistance.begin();imWDistance!=mWDistance.end();imWDistance++)
		vpWDistance.push_back(make_pair((*imWDistance).first, (*imWDistance).second));

	int dThreshold=maxWDistance;

	while(!Heap.empty()){
		Heap.extract_min(topNodeID, topDistance);
		if(topDistance>dThreshold)
			break;

		for(auto ivp=AdjaShort[topNodeID].begin();ivp!=AdjaShort[topNodeID].end();ivp++){
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second;
			if(NodeOrder[neighborNodeID]<=NodeOrder[ID2])
				continue;
			int d = vDistance[topNodeID] + neighborLength;
			if(vDistance[neighborNodeID] ==INF && neighborNodeID != ID1)
			{
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}else if(vDistance[neighborNodeID]>d){
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
			else if(vDistance[neighborNodeID] == d)
			{
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
		}

		for(auto ivp=AdjaShortR[topNodeID].begin();ivp!=AdjaShortR[topNodeID].end();ivp++){
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second;
			if(NodeOrder[neighborNodeID]<=NodeOrder[ID2])
				continue;
			int d = vDistance[topNodeID] + neighborLength;
			if(vDistance[neighborNodeID] ==INF && neighborNodeID != ID1)
			{
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}else if(vDistance[neighborNodeID]>d){
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
			else if(vDistance[neighborNodeID] == d)
			{
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
		}
	}

	for(auto imDistance = mDistance.begin(); imDistance != mDistance.end(); imDistance++)
	{
		if((*imDistance).second > mWDistance[(*imDistance).first])
		{
			int w = (*imDistance).first;
			int distance = mWDistance[(*imDistance).first];

			addedShortcut.push_back(make_pair(w,distance));
		}
	}
	return 0;
}

void Graph::DHLdecStr(int a,int b, int oldW, int newW){

	for(int i=0;i<Neighbor[a].size();i++){
		if(Neighbor[a][i].first==b){
			Neighbor[a][i].second=newW;
			break;
		}
	}
	for(int i=0;i<Neighbor[b].size();i++){
		if(Neighbor[b][i].first==a){
			Neighbor[b][i].second=newW;
			break;
		}
	}

	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC;
	map<pair<int,int>,int> OCdis;

	if(NodeOrder[a]<NodeOrder[b]){
		for(int i=0;i<NeighborCon[a].size();i++){
			if(NeighborCon[a][i].first==b){
				if(NeighborCon[a][i].second.first>newW){
					NeighborCon[a][i].second.first=newW;
					NeighborCon[a][i].second.second=1;

					OCdis[make_pair(a,b)]=newW;
					OC.insert(OrderComp(a,b));
				}else if(NeighborCon[a][i].second.first==newW)
					NeighborCon[a][i].second.second+=1;
				break;
			}
		}
	}else{
		for(int i=0;i<NeighborCon[b].size();i++){
			if(NeighborCon[b][i].first==a){
				if(NeighborCon[b][i].second.first>newW){
					NeighborCon[b][i].second.first=newW;
					NeighborCon[b][i].second.second=1;

					OCdis[make_pair(b,a)]=newW;
					OC.insert(OrderComp(b,a));
				}else if(NeighborCon[b][i].second.first==newW)
					NeighborCon[b][i].second.second+=1;
				break;
			}
		}
	}


	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];
		map<int,int> InM2t; 
		vector<pair<int,int>> InMLower; 
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
				InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
			else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
				InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
		}
		int inID,inW,inWt;
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(InM2t.find(inID)!=InM2t.end()){
				inW=InM2t[inID];
				inWt=NeighborCon[t][i].second.first;
				if(inWt>inW+wt){
					NeighborCon[t][i].second.first=inW+wt;
					NeighborCon[t][i].second.second=1;
					OCdis[make_pair(t,inID)]=inW+wt;
					OrderComp oc={t,inID};
					OC.insert(oc);
				}else if(inWt==inW+wt){
					NeighborCon[t][i].second.second+=1;
				}
			}
		}

		for(int i=0;i<InMLower.size();i++){
			inID=InMLower[i].first; inW=InMLower[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					inWt=NeighborCon[inID][j].second.first;
					if(inWt>inW+wt){
						NeighborCon[inID][j].second.first=inW+wt;
						NeighborCon[inID][j].second.second=1;

						OCdis[make_pair(inID,t)]=inW+wt;
						OrderComp oc={inID,t};
						OC.insert(oc);
					}else if(inWt==inW+wt)
						NeighborCon[inID][j].second.second+=1;
					break;
				}
			}
		}
	}
}

void Graph::DHLincStrMT(int a,int b, int oldW, int newW){

	for(int i=0;i<Neighbor[a].size();i++){
		if(Neighbor[a][i].first==b){
			Neighbor[a][i].second=newW;
			break;
		}
	}
	for(int i=0;i<Neighbor[b].size();i++){
		if(Neighbor[b][i].first==a){
			Neighbor[b][i].second=newW;
			break;
		}
	}

		NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
		set<OrderComp> OC; 
		map<pair<int,int>,int> OCdis;

		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first==oldW){
						NeighborCon[a][i].second.second-=1;
						if(NeighborCon[a][i].second.second<1){
							OrderComp oc={a,b};
							OC.insert(oc);
							OCdis[make_pair(a,b)]=oldW;
						}
					}
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first==oldW){
						NeighborCon[b][i].second.second-=1;
						if(NeighborCon[b][i].second.second<1){
							OrderComp oc={b,a};
							OC.insert(oc);
							OCdis[make_pair(b,a)]=oldW;
						}
					}
					break;
				}
			}
		}


		while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		for(int i=0;i<NeighborCon[s].size();i++){
			inID=NeighborCon[s][i].first;
			inW=NeighborCon[s][i].second.first;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(NeighborCon[t][i].second.first==wt+inW){
					NeighborCon[t][i].second.second-=1;
					if(NeighborCon[t][i].second.second<1){
						OrderComp oc={t,inID};
						OC.insert(oc);
						OCdis[make_pair(t,inID)]=wt+inW;
					}
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					if(NeighborCon[inID][j].second.first==inW+wt){
						NeighborCon[inID][j].second.second-=1;
						if(NeighborCon[inID][j].second.second<1){
							OrderComp oc={inID,t};
							OC.insert(oc);
							OCdis[make_pair(inID,t)]=wt+inW;
						}
					}
					break;
				}
			}
		}

		wt=INF; int countwt=0;
		for(int i=0;i<Neighbor[s].size();i++){
			if(Neighbor[s][i].first==t){
				wt=Neighbor[s][i].second;
				countwt=1;
				break;
			}
		}
		int ssw,wtt,wid;
		vector<int> Wnodes; 
		if(s<t){
			Wnodes=SCconNodesMT[s][t];
		}else{
			Wnodes=SCconNodesMT[t][s];
		}

		for(int i=0;i<Wnodes.size();i++){
			wid=Wnodes[i];
			for(int j=0;j<NeighborCon[wid].size();j++){
				if(NeighborCon[wid][j].first==s){
					ssw=NeighborCon[wid][j].second.first;
				}
				if(NeighborCon[wid][j].first==t){
					wtt=NeighborCon[wid][j].second.first;
				}
			}

			if(ssw+wtt<wt){
				wt=ssw+wtt;
				countwt=1;
			}else if(ssw+wtt==wt){
				countwt+=1;
			}
		}

		for(int i=0;i<NeighborCon[s].size();i++){
			if(NeighborCon[s][i].first==t){
				NeighborCon[s][i].second.first=wt;
				NeighborCon[s][i].second.second=countwt;
				break;
			}
		}
	}
}


void Graph::DHLincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC; 
	map<pair<int,int>,int> OCdis;


	for(int wb=0;wb<wBatch.size();wb++){
		int a=wBatch[wb].first.first;
		int b=wBatch[wb].first.second;
		int oldW=wBatch[wb].second.first;
		int newW=wBatch[wb].second.second;
		for(int i=0;i<Neighbor[a].size();i++){
			if(Neighbor[a][i].first==b){
				Neighbor[a][i].second=newW;
				break;
			}
		}
		for(int i=0;i<Neighbor[b].size();i++){
			if(Neighbor[b][i].first==a){
				Neighbor[b][i].second=newW;
				break;
			}
		}


		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first==oldW){
						NeighborCon[a][i].second.second-=1;
						if(NeighborCon[a][i].second.second<1){
							OrderComp oc={a,b};
							OC.insert(oc);
							OCdis[make_pair(a,b)]=oldW;
						}
					}
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first==oldW){
						NeighborCon[b][i].second.second-=1;
						if(NeighborCon[b][i].second.second<1){
							OrderComp oc={b,a};
							OC.insert(oc);
							OCdis[make_pair(b,a)]=oldW;
						}
					}
					break;
				}
			}
		}
	}

	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());
		wt=OCdis[make_pair(s,t)];
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		for(int i=0;i<NeighborCon[s].size();i++){
			inID=NeighborCon[s][i].first;
			inW=NeighborCon[s][i].second.first;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(NeighborCon[t][i].second.first==wt+inW){
					NeighborCon[t][i].second.second-=1;
					if(NeighborCon[t][i].second.second<1){
						OrderComp oc={t,inID};
						OC.insert(oc);
						OCdis[make_pair(t,inID)]=wt+inW;
					}
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					if(NeighborCon[inID][j].second.first==inW+wt){
						NeighborCon[inID][j].second.second-=1;
						if(NeighborCon[inID][j].second.second<1){
							OrderComp oc={inID,t};
							OC.insert(oc);
							OCdis[make_pair(inID,t)]=wt+inW;
						}
					}
					break;
				}
			}
		}

		wt=INF; int countwt=0;
		for(int i=0;i<Neighbor[s].size();i++){
			if(Neighbor[s][i].first==t){
				wt=Neighbor[s][i].second;
				countwt=1;
				break;
			}
		}
		int ssw,wtt,wid;
		vector<int> Wnodes; 
		if(s<t){
			Wnodes=SCconNodesMT[s][t];
		}else{

			Wnodes=SCconNodesMT[s][t];
		}

		for(int i=0;i<Wnodes.size();i++){
			wid=Wnodes[i];
			for(int j=0;j<NeighborCon[wid].size();j++){
				if(NeighborCon[wid][j].first==s){
					ssw=NeighborCon[wid][j].second.first;
				}
				if(NeighborCon[wid][j].first==t){
					wtt=NeighborCon[wid][j].second.first;
				}
			}

			if(ssw+wtt<wt){
				wt=ssw+wtt;
				countwt=1;
			}else if(ssw+wtt==wt){
				countwt+=1;
			}
		}
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NeighborCon[s][i].first==t){
				NeighborCon[s][i].second.first=wt;
				NeighborCon[s][i].second.second=countwt;
				break;
			}
		}
	}
}


