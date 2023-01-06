#include "head.h"

int main(){
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_cost;
	double runningtime;
	string DesFile="./data2/";
	vector<pair<pair<int,int>,int>> Testdata;
	Graph g;
	string destiGraph=DesFile+"sample";
	string orderfile=DesFile+"dataorder_sample";
	g.ReadGraph(destiGraph);

	t1=std::chrono::high_resolution_clock::now();
	g.DHLconsorder(orderfile);
	t2=std::chrono::high_resolution_clock::now();
	time_cost = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runningtime= time_cost.count();
	cout<<"DHL's label construction time "<<runningtime<<endl;
	g.CorCheckDHL();

	//index size & efficiency
	string ODfile=DesFile+"query";
	g.DHLindexsize();
	g.SCconNodesize();
	g.EffiCheckCH(ODfile);

	//index update data
	string updateFile=DesFile+"update";
	vector<pair<pair<int,int>,int>> testdata;
	g.ReadIns(updateFile, testdata);
	vector<pair<pair<int,int>,pair<int,int>>> testdataInc, testdataDec;
	for(int k=0;k<testdata.size();k++){
		testdataInc.push_back(make_pair(testdata[k].first,make_pair(testdata[k].second,testdata[k].second*0.8)));
		testdataDec.push_back(make_pair(testdata[k].first,make_pair(testdata[k].second,testdata[k].second*0.8)));
	}


	//DHL for weight update
	Graph gu=g;
	t1=std::chrono::high_resolution_clock::now();
	gu.DHLincBatMT(testdataInc);
	t2=std::chrono::high_resolution_clock::now();
	time_cost = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runningtime= time_cost.count();
	cout<<"query for weight update(increase) "<<runningtime<<" "<<runningtime/testdata.size()<<endl;

	return 0;
}



