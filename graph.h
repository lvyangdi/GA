/*
 * Authors: Yangdi Lyu
 * Jan 2018
 */

#ifndef __GRAPH_H__
#define __GRAPH_H__
#include <string>
#include <vector>
using namespace std;

class Graph;
class Edge;
class Vertex;

class Graph 
{
public:
    string          header;
    vector<Vertex*>  vertexArr;      /* Store all vertices */
    vector<Edge*>    edgArr;         /* Store all edges */
    vector<Vertex*> combGates;
    vector<Vertex*> seqGates;
    vector<Edge*>   primInEdges;   /* Vector of Primary Inputs */
    vector<Edge*>   primOutEdges;  /* Vector of Primary Outputs */
    vector<Edge*>   lowprobEdges;  /* Vector of rare edges */
    vector<Vertex*> topoArr;

    Graph(string name = "") {header = name;}
    ~Graph();

    Edge* findEdge(const string& name);
    bool getTopoOrder();
    Edge* createEdge(string& name);
    Vertex* createNewVertexForVlog(vector<string>& tokens);
    void populateOnlyNewValues(const string& vec);    // Only newVal will be changed
    int populatePairValues(const string& vec1, const string& vec2);       // Both oldVal and newVal will be changed
    int populateValues(const string& vec);           // Both oldVal and newVal will be changed
    int countTransitions();                     // Compare old values with new values, compute the diff
    
    void populateAfterFlipping(Edge* edg);
};

class Vertex
{
public:
    int idx;                 /* Assigned when topo sorting */
    int bPrimIn;             /* if it's a primary input or not */
    int bPrimOut;
    string name;
    string type;
    bool isSeq;
    int level;      /* It specifies the level of gate */  
    vector<Edge*> inEdges;
    vector<Edge*> outEdges;
    vector<string> inPorts;
    vector<string> outPorts;

    float  nInput00;
    float  nInput01;
    float  nInput10;
    float  nInput11;

    Vertex(string tp, string nm):
        bPrimIn(0), bPrimOut(0), name(nm), type(tp), isSeq(false), 
        level(0), nInput00(0), nInput01(0), nInput10(0),
        nInput11(0){}

    Edge* findDataEdge();
    int calculateOutTransition();
    int calculateOutTransitionNewOnly();
    static void getValue(string& type, vector<int>& DIN, int& val1, int& val2);
};

class Edge
{
public:
    string name;
    int oldVal;
    int newVal;
    int lowProbVal;
    bool bPrimIn;
    bool bPrimOut;
    Vertex* fromNode;
    vector<Vertex*> toNodes;

    float    signal0Prob;
    float    signal1Prob;
    float    activity01;
    float    activity10;

    Edge(string nm):
        name(nm), oldVal(0), newVal(1), lowProbVal(0), bPrimIn(false),
        bPrimOut(false), fromNode(nullptr), signal0Prob(0), 
        signal1Prob(0), activity01(0), activity10(0){}
};

Graph* createGraphFromVerilog(string filename, string headername, int& primInCnt);

void tokenize(string line, vector<string>& tokens);

#endif
