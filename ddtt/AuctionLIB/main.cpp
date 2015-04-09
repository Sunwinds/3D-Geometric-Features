
///////////////////////////////////////////////////////////////////////////////
// File name: main.cpp
// First obtain assignment, then call the algorithm.
// Lantao Liu, Apr 15, 2013
/////////////////////////////////////////////////////////////////////////////// 

#include "Define.h"
#include "CmdParser.h"
#include "Assignment.h"
#include "iAuction.h"
#include "Utils.h"

#include "Hungarian.h"
#include "BipartiteGraph.h"

using namespace AssignmentLib;

int hungarian_main(int argc, char** argv)
{
    //define a matrix;
    Matrix m;

    //parse command line and generate an assignment
    CmdParser parser(argc, argv);
    parser.DisplayInput();

    Assignment as;

    /*if(parser.GetInputFileName().empty()){
        if(parser.GetSeed())
            as.SetSeed(parser.GetSeed());
        //else
        //    as.SetSeed(time(NULL));
        cout<<endl<<"  *Seed for random generator: "<<as.GetSeed()<<endl;
        m=as.RandomGenerate( parser.GetAssignmentSize(),
                             parser.GetAssignmentSize(),
                             MAX_RANDOM, 0 );
    }
    else{
        ifstream myfile(parser.GetInputFileName().c_str());
        m=as.ImportMatrix(myfile);

    }
    as.DisplayMatrix(m);*/

    //define a bipartite graph
    BipartiteGraph bg(m);

    //run Hungarian method
    Hungarian h(bg);
    h.HungarianAlgo();

    //for testing
    //Test t;
    //t.testHungarian();

    return 0;
}


int library_main(int argc, char** argv)
{

    //parse command line and generate an assignment
    CmdParser parser(argc, argv);
    parser.DisplayInput();

    Assignment as;
    if(parser.GetInputFileName().empty()){
        if(parser.GetSeed())
            as.SetSeed(parser.GetSeed());
#ifndef WIN32
        else
            as.SetSeed(time(NULL));
#endif

        _cout(endl<<"  *Seed for random generator: "<<as.GetSeed()<<endl);
        _cout("  RAND_MAX: "<<RAND_MAX<<endl);
        as.RandomGenerate( parser.GetAssignmentSize(),
                           parser.GetAssignmentSize(),
                           MAX_RANDOM,
                           parser.GetSeed() );
    }
    else{
        ifstream myfile(parser.GetInputFileName().c_str());
        as.ImportMatrix(myfile);
    }
    as.DisplayMatrix(as.GetMatrix());

    Utils utils;

    try{

        iAuction auction(as.GetMatrix());

        double time_start = utils.GetCurTime();
        double time = 0;

        auction.MainAlgo();

        time += utils.GetCurTime() - time_start;
        cout<<"Time used: "<<time<<endl;
    }
    catch(int e){

        if(e==EXCEPTION_BROKEN){
            cerr<<"I suspect there is something wrong..."<<endl;
            cerr<<"Current seed: "<<as.GetSeed()<<endl;
            cerr<<"The problem has been written in to \"debug.txt\"."<<endl;
            utils.WriteMatrix(as);
        }
        else if(e==EXCEPTION_WRONG){
            cerr<<"!!Wrong solution, double check it. "<<endl;
            cerr<<"Current seed: "<<as.GetSeed()<<endl;
            cerr<<"The problem has been written in to \"debug.txt\"."<<endl;
            utils.WriteMatrix(as);
        }
        else
            cerr<<"Unknown exception caught. "<<endl;
    }

    return 0;

}


