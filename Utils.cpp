#include "head.h"

void Graph::NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
	sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		for(int h=0;h<Neighvec.size();h++){
			ID2=Neighvec[h].first;
			insertEMTOrderGenerate(ID1, ID2, 1);
		}
	}
	sm->notify();
}

void Graph::insertEMTOrderGenerate(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
		DD[u]++;
		DD2[u]++;
	}
}

void Graph::NoRegularGraphOrder(){
	cout<<"begin order the degree!"<<endl;
	benchmark::heap<2, int, int> Q(nodenum);
	for(int i=0;i<nodenum;i++){
		Q.update(i,Neighbor[i].size());
		if(Neighbor.size()==0) cout<<"Isolated node "<<i<<endl;
	}
	int topID, topdegree;
	int cnt=0;
	NodeOrder.assign(nodenum, 0);
	vNodeOrder.assign(nodenum, 0);
	while(!Q.empty()){
		Q.extract_min(topID, topdegree);
		NodeOrder[topID]=cnt;
		vNodeOrder[cnt]=topID;
		cnt+=1;
	}
	cout<<"Node finish ordering!"<<endl;
}

void Graph::ReadGraph(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}

	int nn,en;
	IF>>nodenum>>edgenum;

	eNum=0;
	set<pair<int,int>> eSet;

	vector<pair<int,int>> vecp; vecp.clear();
	Neighbor.assign(nodenum, vecp);

	set<int> setp; setp.clear();
	AdjacentNodes.assign(nodenum, setp);
	set<pair<int,int>> EdgeRedun;

	int ID1, ID2, weight;
	for(int i=0;i<edgenum;i++){
		IF>>ID1>>ID2>>weight;
		if(EdgeRedun.find(make_pair(ID1,ID2))==EdgeRedun.end()){
			Neighbor[ID1].push_back(make_pair(ID2, weight));
			AdjacentNodes[ID1].insert(ID2);
		}
		EdgeRedun.insert(make_pair(ID1,ID2));
	}
}

void Graph::EffiCheckDij(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}


	int s, t, d;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<200;i++){
		s=ODpair[i].first;
		t=ODpair[i].second;
		d=BiDij(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"BiDij query time "<<runT/ODpair.size()<<endl;
}


void Graph::EffiCheckCH(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}

	int s, t, d;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair.size();i++){
		s=ODpair[i].first;
		t=ODpair[i].second;
		d=QueryDHL(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"query time (both border)"<<runT/ODpair.size()<<endl;
}


void Graph::CorCheckDij(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dij(s,t);
		d2=BiDij(s,t);
		if(d1!=d2)
		    cout<<"Correct"<<endl;
	}
}

void Graph::CorCheckDHL(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dij(s,t);
		d2=QueryDHL(s,t);
		if(d1!=d2)
		   cout<<"Build TDIF between"<<" "<<s<<" "<<t<<endl;
	}
}

int Graph::Dij(int ID1, int ID2){
	if(ID1==ID2) return 0;
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);

	distance[ID1]=0;
	int topNodeID, topNodeDis;
	int NNodeID,NWeigh;

	int d=INF;

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		if(topNodeID==ID2){
			d=distance[ID2];
			break;
		}
		closed[topNodeID]=true;

		for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
				}
			}
		}
	}
	return d;
}

int Graph::BiDij(int ID1, int ID2){
	if(ID1==ID2) return 0;
	benchmark::heap<2, int, int> queueF(nodenum), queueB(nodenum);
	queueF.update(ID1,0);
	queueB.update(ID2,0);

	vector<bool> closedF(nodenum, false), closedB(nodenum, false);
	vector<int> distanceF(nodenum, INF), distanceB(nodenum, INF);

	distanceF[ID1]=0;
	distanceB[ID2]=0;
	int topNodeIDF, topNodeDisF, topNodeIDB, topNodeDisB;
	int NNodeIDF,NWeighF, NNodeIDB, NWeighB;
	int d=INF;

	while(!queueF.empty() || !queueB.empty()){
		if(queueF.top()+queueB.top()>=d){
			return d;
		}
		queueF.extract_min(topNodeIDF, topNodeDisF);
		closedF[topNodeIDF]=true;
		for(auto itF=Neighbor[topNodeIDF].begin();itF!=Neighbor[topNodeIDF].end();itF++){
			NNodeIDF=(*itF).first;
			NWeighF=(*itF).second+topNodeDisF;
			if(closedB[NNodeIDF] && NWeighF+distanceB[NNodeIDF]<d){
				d=NWeighF+distanceB[NNodeIDF];
			}
			if(!closedF[NNodeIDF]){
				if(distanceF[NNodeIDF]>NWeighF){
					distanceF[NNodeIDF]=NWeighF;
					queueF.update(NNodeIDF, NWeighF);
				}
			}
		}
		queueB.extract_min(topNodeIDB, topNodeDisB);
		closedB[topNodeIDB]=true;
		for(auto itB=Neighbor[topNodeIDB].begin();itB!=Neighbor[topNodeIDB].end();itB++){
			NNodeIDB=(*itB).first;
			NWeighB=(*itB).second+topNodeDisB;
			if(closedF[NNodeIDB] && NWeighB+distanceF[NNodeIDB]<d){
				d=NWeighB+distanceF[NNodeIDB];
			}
			if(!closedB[NNodeIDB]){
				if(distanceB[NNodeIDB]>NWeighB){
					distanceB[NNodeIDB]=NWeighB;
					queueB.update(NNodeIDB, NWeighB);
				}
			}
		}
	}
	return d;
}


int Graph::Astar(int ID1, int ID2){
	if(ID1==ID2) return 0;
	benchmark::heap<2, int, int> pqueue(nodenum);
	int heurisDis=EuclideanDis(ID1,ID2);
	pqueue.update(ID1,heurisDis);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);

	distance[ID1]=0;
	int topNodeID, topNodeDis;
	int NNodeID,NWeigh;
	int d=INF;
	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		if(topNodeID==ID2){
			d=distance[ID2];
			break;
		}
		closed[topNodeID]=true;

		for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+distance[topNodeID];
			heurisDis=EuclideanDis(NNodeID,ID2);
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh+heurisDis);
				}
			}
		}
	}
	return d;
}

int Graph::EuclideanDis(int s, int t){
	int lat=(int)(abs(GraphLocation[s].first-GraphLocation[t].first)*111319);
	int lon=(int)(abs(GraphLocation[s].second-GraphLocation[t].second)*83907);
	int min,max;
	min=(lat>lon)?lon:lat;
	max=(lat>lon)?lat:lon;
	int approx=max*1007+min*441;
	if(max<(min<<4))
		approx-=max*40;
	return (approx+512)>>10;
}

int Graph::writeShortCutorder(string filename){
	ofstream ofile(filename);
	cout<<"Writing shortcuts!"<<endl;
	for(int i = 0; i < nodenum; i++)
	{
		ofile << i << "\t" << NodeOrder[i] << "\t";

		vector<pair<int, int>> shortcut;
		int ID, Wei;
		for(int k=0;k<vvpShortCut[i].size();k++){
			ID=vvpShortCut[i][k].first;
			Wei=vvpShortCut[i][k].second;
			if(NodeOrder[ID]>NodeOrder[i])
				shortcut.push_back(make_pair(ID, Wei));
		}

		for(int j=0;j<Neighbor[i].size();j++){
			ID=Neighbor[i][j].first;
			Wei=Neighbor[i][j].second;
			if(NodeOrder[ID]>NodeOrder[i])
				shortcut.push_back(make_pair(ID, Wei));
		}

		if(shortcut.size() != 0)
		{
			ofile << 1;
			ofile << "\t" << shortcut.size();
			vector<pair<int,int>>::iterator iR;
			iR=shortcut.begin();
			while(iR!=shortcut.end()){
				ofile<<"\t"<<(*iR).first<<"\t"<<(*iR).second;
				iR++;
			}
		}
		else
			ofile << 0;

		ofile << endl;
	}
	return 0;
}

int Graph::ReadShortCut(string filename){
	ifstream inSCH(filename);
	int nodeID, nodeOrder, c, nsc, nrsc, n, d, m;
	cout<<"Reading shortcuts!"<<endl;

	for(int i = 0; i < nodenum; i++)
	{
		inSCH >> nodeID >> nodeOrder >> c;
		NodeOrder[nodeID] = nodeOrder;
		if(!c)
			continue;

		inSCH >> nrsc;
		for(int j = 0; j < nrsc; j++)
		{
			inSCH >> n >> d;
			AdjaShort[nodeID].push_back(make_pair(n,d));
			AdjaShortR[n].push_back(make_pair(nodeID,d));
		}
	}
	inSCH.close();

	cout<<"shortcut finish reading!"<<endl;
	return 0;
}

void Graph::writeDHLIncrease(string indexfile)
{
	ofstream ofSupportNodes(indexfile+"SupportNodes");
	ofstream ofPathInfor(indexfile+"PathInfor");
	ofstream ofEdgeOnPath(indexfile+"EdgeOnPath");

	for(int i = 0; i < nodenum; i++)
	{
		ofSupportNodes << i << "\t" << SupportNodes[i].size() << endl;
		for(auto& im : SupportNodes[i])
		{
			ofSupportNodes << im.first << "\t" << im.second.size();
			for(auto& iv : im.second)
				ofSupportNodes << "\t" << iv;
			ofSupportNodes << endl;
		}
	}
	ofSupportNodes.close();

	for(int i = 0; i < nodenum; i++)
	{
		ofPathInfor << i << "\t" << PathInfor[i].size() << endl;
		for(auto& im : PathInfor[i])
			ofPathInfor << im.first.first << "\t" << im.first.second << "\t" << im.second << endl;
	}
	ofPathInfor.close();

	for(int i = 0; i < (int)EdgeOnPath.size(); i++)
	{
		ofEdgeOnPath << i << "\t" << EdgeOnPath[i].size() << endl;
		for(auto& iTri : EdgeOnPath[i])
			ofEdgeOnPath << iTri.u << "\t" << iTri.v << "\t" << iTri.w << endl;
	}
	ofEdgeOnPath.close();
}

void Graph::readDHLIncrease(string indexfile)
{
	set<pair<int,int>> setpair;
	InvalidWP.assign(nodenum, setpair);
	ifstream ifSupportNodes(indexfile + "SupportNodes");
	ifstream ifPathInfor(indexfile + "PathInfor");
	ifstream ifEdgeOnPath(indexfile + "EdgeOnPath");

	int nID, ssize;
	for(int i = 0; i < nodenum; i++)
	{
		ifSupportNodes >> nID  >> ssize;
		int n, nsize;
		map<int, vector<int> > mv;
		for(int j = 0; j < ssize; j++)
		{
			ifSupportNodes >> n >> nsize;
			vector<int> vSupport(nsize, -1);
			for(int k = 0; k < nsize; k++)
				ifSupportNodes >> vSupport[k];
			mv[n] = vSupport;
		}

		SupportNodes[i] = mv;
	}
	ifSupportNodes.close();

	for(int i = 0; i < nodenum; i++)
	{
		ifPathInfor >> nID >> ssize;
		unordered_map<pair<int,int>, int, hash_pair> pathinfvec;
		int a, b, c;
		for(int j = 0; j < ssize; j++)
		{
			ifPathInfor >> a >> b >> c;
			pathinfvec.insert(make_pair(make_pair(a, b), c));
		}
		PathInfor[i] = pathinfvec;
	}
	ifPathInfor.close();

	int eID;
	for(int i = 0; i < (int)EdgeOnPath.size(); i++)
	{
		ifEdgeOnPath >> eID >> ssize;
		tri t;
		vector<tri> vt(ssize, t);
		for(int j = 0; j < ssize; j++)
		{
			tri t2;
			ifEdgeOnPath >> t2.u >> t2.v >> t2.w;
			vt[j] = t2;
		}
		EdgeOnPath[i] = vt;
	}
	ifEdgeOnPath.close();
}

int Graph::writeShortCutDouble(string filename){

	ofstream ofile(filename);
	cout<<"Writing shortcuts!"<<endl;
	for(int i = 0; i < nodenum; i++)
	{
		ofile << i << "\t" << NodeOrder[i] << "\t";

		if(NeighborCon[i].size() != 0)
		{
			ofile << 1;

			ofile << "\t" << NeighborCon[i].size();
			vector<pair<int,pair<int,int>>>::iterator iR;
			iR=NeighborCon[i].begin();
			while(iR!=NeighborCon[i].end()){
				ofile<<"\t"<<(*iR).first<<"\t"<<(*iR).second.first<<"\t"<<(*iR).second.second;
				iR++;
			}
		}
		else
			ofile << 0;

		ofile << endl;
	}
	return 0;
}

int Graph::ReadShortCutDouble(string filename){
	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);
	ifstream inSCH(filename);
	int nodeID, nodeOrder, c, nsc, nrsc, n, d, m;
	for(int i = 0; i < nodenum; i++)
	{
		inSCH >> nodeID >> nodeOrder >> c;
		NodeOrder[nodeID] = nodeOrder;
		if(!c)
			continue;

		inSCH >> nrsc;
		for(int j = 0; j < nrsc; j++)
		{
			inSCH >> n >> d >> m;
			NeighborCon[nodeID].push_back(make_pair(n,make_pair(d,m)));
		}
	}
	inSCH.close();

	cout<<"Finish reading!"<<endl;
	return 0;
}

void Graph::ReadIns(string filename,vector<pair<pair<int,int>,int>>& TestData){
	TestData.clear();

	int num, ID1, ID2, neww;
	ifstream IF(filename);
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>ID1>>ID2>>neww;
		TestData.push_back(make_pair(make_pair(ID1, ID2), neww));
	}
	IF.close();
}


long long Graph::DHLindexsize(){
	long long m=0;

	for(int i=0;i<NeighborCon.size();i++){
		m+=NeighborCon[i].size()*3*sizeof(int);
	}

	cout<<"DHL index size(border labeling list) "<<(double)m/1024/1024<<endl;
	return m;
}


long long Graph::SCconNodesize(){
	long long m=0;

	for(int i=0;i< SCconNodesMT.size();i++){
		for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
			m+=sizeof(int)+(*it).second.size()*sizeof(int);
		}
	}

	cout<<"DHL Support Nodes size "<<(double)m/1024/1024<<endl;
	return m;
}



int	Graph::QueryDHLorder(int ID1, int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

	int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

	while(!fHeapForward.empty() || !fHeapBackward.empty() )
	{
		if(bF && bB)
			break;
		if(bF && fHeapBackward.empty())
			break;
		if(bB && fHeapForward.empty())
			break;
		if(!fHeapForward.empty() && !bF)
		{
			fHeapForward.extract_min(topNodeIDForward, topDisForward);

			if(vDistanceForward[topNodeIDForward] > d)
				bF = true;

			vVisitedF[topNodeIDForward] = true;

			if(vVisitedB[topNodeIDForward]){
				int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
				if(distTmp<d){
					d=distTmp;
				}
			}

			for(auto out=AdjaShort[topNodeIDForward].begin();out!=AdjaShort[topNodeIDForward].end();out++){
				neighborNodeID = (*out).first;
				neighborLength = (*out).second;

				int df = vDistanceForward[topNodeIDForward] + neighborLength;
				if(!vVisitedF[neighborNodeID]){
					if(vDistanceForward[neighborNodeID] > df){
						vDistanceForward[neighborNodeID] = df;
						fHeapForward.update(neighborNodeID, df);
					}
				}
			}
		}

		if(!fHeapBackward.empty() && !bB)
		{
			fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

			if(vDistanceBackward[topNodeIDBackward] > d)
				bB = true;

			vVisitedB[topNodeIDBackward] = true;

			if(vVisitedF[topNodeIDBackward]){
				int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
				if(distTmp<d){
					d=distTmp;
				}
			}

			for(auto in=AdjaShort[topNodeIDBackward].begin();in!=AdjaShort[topNodeIDBackward].end();in++){
				neighborNodeID = (*in).first;
				neighborLength = (*in).second;

				int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
				if(!vVisitedB[neighborNodeID]){
					if(vDistanceBackward[neighborNodeID]>db){
						vDistanceBackward[neighborNodeID] = db;
						fHeapBackward.update(neighborNodeID, db);
					}
				}
			}
		}
	}
	return d;
}

int	Graph::QueryDHL(int ID1, int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

	int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

	while(!fHeapForward.empty() || !fHeapBackward.empty() )
	{
		if(bF && bB)
			break;
		if(bF && fHeapBackward.empty())
			break;
		if(bB && fHeapForward.empty())
			break;
		if(!fHeapForward.empty() && !bF)
		{
			fHeapForward.extract_min(topNodeIDForward, topDisForward);

			if(vDistanceForward[topNodeIDForward] > d)
				bF = true;

			vVisitedF[topNodeIDForward] = true;

			if(vVisitedB[topNodeIDForward]){
				int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
				if(distTmp<d){
					d=distTmp;
				}
			}

			for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
				neighborNodeID = (*out).first;
				neighborLength = (*out).second.first;

				int df = vDistanceForward[topNodeIDForward] + neighborLength;
				if(!vVisitedF[neighborNodeID]){
					if(vDistanceForward[neighborNodeID] > df){
						vDistanceForward[neighborNodeID] = df;
						fHeapForward.update(neighborNodeID, df);
					}
				}
			}
		}


		if(!fHeapBackward.empty() && !bB)
		{
			fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

			if(vDistanceBackward[topNodeIDBackward] > d)
				bB = true;

			vVisitedB[topNodeIDBackward] = true;

			if(vVisitedF[topNodeIDBackward]){
				int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
				if(distTmp<d){
					d=distTmp;
				}
			}

			for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
				neighborNodeID = (*in).first;
				neighborLength = (*in).second.first;

				int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
				if(!vVisitedB[neighborNodeID]){
					if(vDistanceBackward[neighborNodeID]>db){
						vDistanceBackward[neighborNodeID] = db;
						fHeapBackward.update(neighborNodeID, db);
					}
				}
			}
		}
	}
	return d;
}

int Graph::ShortestDisQuery(int ID1,int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(Label[ID2].find(hub)!=Label[ID2].end()){
			dis2=Label[ID2][hub];
			if(dis1+dis2<d){
				d=dis1+dis2;
			}
		}
	}
	return d;
}
