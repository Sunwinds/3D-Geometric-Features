
#include "PlotGraph.h"
using namespace AssignmentLib;

#include <sstream>
#include <iomanip>

//#define SAVE_PLOTS
int Cnt = 0;

void 
PlotGraph::PlotBipartiteGraph(BipartiteGraph& _bg,
			      vector<VID>& _S, 
			      vector<VID>& _T,
			      vector<VID>& _N,
			      vector<EID>& _EG,
			      vector<EID>& _M, 
			      int target_task){


}


void
PlotGraph::PlotAugmentingPath(BipartiteGraph& _bg, vector<EID>& _path){

}

void
PlotGraph::DisplayData(const vector<VID>& vs){

  for(vector<VID>::const_iterator itr = vs.begin(); itr != vs.end(); itr++)
    cout<<*itr<<" ";
  cout<<endl;

}

void
PlotGraph::DisplayData(const vector<EID>& es){

  for(vector<EID>::const_iterator itr = es.begin(); itr != es.end(); itr++)
    cout<<"("<<itr->first<<","<<itr->second<<") ";
  cout<<endl;

}







