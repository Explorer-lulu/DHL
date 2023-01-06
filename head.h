#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <boost/thread/thread.hpp>
#include "Flag.h"

#include <functional>
#include <bits/stdc++.h>
#include <utility>

struct hash_pair {
   template <class T1, class T2>
   size_t operator()(const pair<T1, T2>& p) const{
      auto hash1 = hash<T1>{}(p.first);
      auto hash2 = hash<T2>{}(p.second);
      return hash1 ^ hash2;
   }
};

#define INF 999999999
using namespace std;
using namespace boost;

struct Nei{            //Neighbor
	int nid;
	int w;
	int c;
};

struct tri{
	int u;
	int v;
	int w;
};

struct Node{
	vector<pair<int,pair<int,int>>> vert;
	vector<pair<int,Nei>> neighInf;
	vector<int> pos, pos2;
	vector<int> dis, cnt;
	set<int> changedPos;
	vector<bool> FN;
	set<int> DisRe;
	vector<int> ch;   
	int height, hdepth;
	int pa;
	int uniqueVertex;
	vector<int> piv;
	Node(){
		vert.clear();
		neighInf.clear();
		pos.clear();
		dis.clear();
		cnt.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
		changedPos.clear();
		FN.clear();
		DisRe.clear();
		piv.clear();
	}
};

class Graph{
public:
	int nodenum;
	int edgenum;

	int threadnum=20;

	int eNum;
	vector<pair<int,int>> Edge;
	unordered_map<pair<int,int>, int, hash_pair> EdgeRe;
	vector<set<int>> E1;


	vector<vector<pair<int,int>>> Neighbor;
	vector<pair<double,double>> GraphLocation;


	vector<int> NodeOrder;
	vector<int> vNodeOrder;

	void HGPbuild();
	void DHLGenerate();
	void hierarchy_shortestd_for_node_load();
	void insertEMTOrderGenerate(int u,int v,int w);
    void topdownLabelupdate();
	void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
	void NoRegularGraphOrder();
	void ReadGraph(string filename);

	void ReadIns(string filename, vector<pair<pair<int,int>,int>>& data);
	void CorCheckDij();
	void CorCheckDHL();
	void EffiCheckDij(string filename);
	void EffiCheckCH(string filename);

	int	QueryDHLorder(int ID1, int ID2);
	int	QueryDHL(int ID1, int ID2);
	int ShortestDisQuery(int ID1,int ID2);


	int Dij(int ID1, int ID2);
	int BiDij(int ID1, int ID2);
	int Astar(int ID1, int ID2);
	int EuclideanDis(int s, int t);


	long long DHLindexsize();
	long long SCconNodesize();


	int DHLcontractionlevel(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult);

	int StartCHPOrderMTWP(string indexfile, string orderfile);
	void CHPconsorderMTWP();
	int DHLsyncWP(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult);

	int writeShortCutorder(string filename);
	int ReadShortCut(string filename);
	void writeDHLIncrease(string indexfile);
	void readDHLIncrease(string indexfile);


	Flag* smLock = new Flag(1);
	vector<map<int,int>> OutEdgesM;
	vector<vector<pair<int, int>>> vvNode;
	vector<vector<pair<int, int>>> vvpShortCut;
	vector<vector<pair<int,pair<int,int>>>> OutEdgesV;
	vector<int>	vnContractedNeighbors;
	vector<int> vnContractedNeighborsOld;
	vector<vector<pair<int, int>>> AdjaShort;
	vector<map<int, vector<int>>> SupportNodes;
	vector<unordered_map<pair<int,int>, int, hash_pair>> PathInfor;
	vector<vector<tri>> EdgeOnPath;
	vector<vector<pair<int, int>>> AdjaShortR;
	vector<Flag*> vSm;
	Flag* sm = new Flag(threadnum);
	int DHLcontractInc(int ID1, int ID2, int dUV, vector<pair<int, int> >& vW, vector<pair<int,int>>& addedShortcut);
	int NewSCweight(int s, int t);
	vector<set<pair<int,int>>> InvalidWP;


	int StartDHL(string indexfile, string orderfile);
	void DHLconsorder(string orderfile);
	void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
	void insertEdge(int u,int v,int w);

	int writeShortCutDouble(string filename);
	int ReadShortCutDouble(string filename);
	//index update
	void DHLdecStr(int a,int b, int oldW, int newW);
	void DHLincStrMT(int a,int b, int oldW, int newW);
	void DHLincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//DHL update


	vector<int> DD,DD2;
	vector<map<int,pair<int,int>>> E;
	vector<vector<pair<int,pair<int,int>>>> NeighborCon;
	vector<map<int, vector<int>>> SCconNodesMT;
	void deleteEorder(int u,int v);
	void insertEorder(int u,int v,int w);
	void deleteE(int u,int v);


	//index update


	vector<int> rank;
	vector<Node> Tree;
	int heightMax;

	vector<vector<int>> VidtoTNid;
	vector<int> EulerSeq;
	vector<int> toRMQ;
	vector<vector<int>> RMQIndex;



	vector<unordered_map<int,int>> Label;
	vector<unordered_map<int,vector<int>>> PruningPointNew;
	set<pair<int,int>> NoSupportedPair;
	int DisQueryPeak(int ID1, int ID2);
	int DisQueryVally(int ID1, int ID2);
	int DisQueryLower1(int ID1, int ID2);
	int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
	void DijkstraPrune1(int NodeID, int ID, int d);
	int PrefixalDisQuery(int ID1, int ID2);
	int DijkstraList(int ID1, vector<int>& distance);

	pair<bool, pair<int,int>> VertexPair(int ID, int hopnum);


	vector<int> BiDijPath(int ID1, int ID2);
	vector<int> AstarPath(int ID1, int ID2);

	void CHconsPath(string orderfile);
	void NeighborComPath(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
	bool insertEPath(int u,int v,int w);
	bool insertEMTPath(int u,int v,int w);

	void makeIndexPath();
	void makeIndexDFSPath(int p, vector<int>& list);

	vector<set<int>> AdjacentNodes;
	vector<vector<pair<int,int>>> label;
	vector<vector<pair<int,pair<int,int>>>> labelPre;
	vector<unordered_map<int,int>> labelPos;
	int DistanceCompute(int s, int t);
	vector<int> PointCompute(int ID1, int ID2);


};

namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

template<int log_k, typename k_t, typename id_t>
class heap {

public:

	typedef k_t key_t;
	typedef id_t node_t;

	static const node_t k = 1 << log_k;

	struct element_t {
		key_t key;
		node_t element;
		element_t() : key(0), element(0) {}
		element_t(const key_t k, const node_t e) : key(k), element(e) {}
	};


public:


	heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {
	}

	heap() {

	}

	inline node_t size() const {
		return n;
	}


	inline bool empty() const {
		return size() == 0;
	}


	inline void extract_min(node_t &element, key_t &key) {
		assert(!empty());

		element_t &front = elements[0];


		element = front.element;
		key = front.key;

		position[element] = NULLINDEX;
		--n;
		if (!empty()) {
			front = elements[n];
			position[front.element] = 0;
			sift_down(0);
		}
	}

	inline key_t top() {
		assert(!empty());

		element_t &front = elements[0];

		return front.key;

	}

	inline node_t top_value() {

		assert(!empty());

		element_t &front = elements[0];

		return front.element;
	}


	inline void update(const node_t element, const key_t key) {
		if (position[element] == NULLINDEX) {
			element_t &back = elements[n];
			back.key = key;
			back.element = element;
			position[element] = n;
			sift_up(n++);
		} else {
			node_t el_pos = position[element];
			element_t &el = elements[el_pos];
			if (key > el.key) {
				el.key = key;
				sift_down(el_pos);
			} else {
				el.key = key;
				sift_up(el_pos);
			}
		}
	}



	inline void clear() {
		for (node_t i = 0; i < n; ++i) {
			position[elements[i].element] = NULLINDEX;
		}
		n = 0;
	}


	inline void clear(node_t v) {
		position[v] = NULLINDEX;
	}

	inline void clear_n() {
		n = 0;
	}



	inline bool contains(const node_t element) const {
		return position[element] != NULLINDEX;
	}


protected:


	inline void sift_up(node_t i) {
		assert(i < n);
		node_t cur_i = i;
		while (cur_i > 0) {
			node_t parent_i = (cur_i-1) >> log_k;
			if (elements[parent_i].key > elements[cur_i].key)
				swap(cur_i, parent_i);
			else
				break;
			cur_i = parent_i;
		}
	}

	inline void sift_down(node_t i) {
		assert(i < n);

		while (true) {
			node_t min_ind = i;
			key_t min_key = elements[i].key;

			node_t child_ind_l = (i << log_k) + 1;
			node_t child_ind_u = std::min(child_ind_l + k, n);

			for (node_t j = child_ind_l; j < child_ind_u; ++j) {
				if (elements[j].key < min_key) {
					min_ind = j;
					min_key = elements[j].key;
				}
			}

			if (min_ind != i) {
				swap(i, min_ind);
				i = min_ind;
			} else {
				break;
			}
		}
	}


	inline void swap(const node_t i, const node_t j) {
		element_t &el_i = elements[i];
		element_t &el_j = elements[j];


		position[el_i.element] = j;
		position[el_j.element] = i;


		element_t temp = el_i;
		el_i = el_j;
		el_j = temp;
	}



private:


	node_t n;

	node_t max_n;

	vector<element_t> elements;

	vector<node_t> position;
};
}

#endif /* HEAD_H_ */
