/*
 * Authors: Yangdi Lyu
 * Jan 2018
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstddef>        // std::size_t
#include <sstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <string.h>
#include "graph.h"

Graph::~Graph() {
    for (Vertex *v : vertexArr)
        delete v;
    for (Edge* e : edgArr)
        delete e;
}

/*
 * Use vec to simulate the graph
 * Only newVal will be changed
 * OldVal is not changed
*/
void Graph::populateOnlyNewValues(const string& vec) {
    size_t noscansize = primInEdges.size();
    // Full scan
    if (vec.length() > noscansize) {
        // Flipflop get value from vec
        for(size_t i = 0; i < seqGates.size(); i++)
        {
            Vertex* s = seqGates[i];
            for (size_t j = 0; j < s->outEdges.size(); j++)
            {
                Edge* edg2 = s->outEdges[j];
                string port = s->outPorts[j];
                if (port.compare(".Q") == 0) {
                    edg2->newVal = (vec[noscansize + i] == '0' ? 0 : 1);
                }
                else if (port.compare(".QN") == 0) {
                    edg2->newVal = (vec[noscansize + i] == '0' ? 1 : 0);
                }
            }
        }
    }
    // no scan
    else {
        //
        for(size_t i = 0; i < seqGates.size(); i++)
        {
            Vertex* s = seqGates[i];
            Edge* edg1 = s->findDataEdge();
            for (size_t j = 0; j < s->outEdges.size(); j++)
            {
                Edge* edg2 = s->outEdges[j];
                string port = s->outPorts[j];
                if (port.compare(".Q") == 0) {
                    edg2->newVal = edg1->newVal;
                }
                else if (port.compare(".QN") == 0) {
                    edg2->newVal = (edg1->newVal == 0 ? 1 : 0);
                }
            }
        }
    }

    for (size_t i = 0; i < noscansize; i++) {
        primInEdges[i]->newVal = (vec[i] == '0') ? 0 : 1;
    }

    for (Vertex* v : topoArr)
        v->calculateOutTransitionNewOnly();
}

/*
 * Use vec1 and vec2 to simulate the graph
 * both oldVal newVal will be changed
 * oldVal from vec1 and newVal from vec2
*/ 
int Graph::populatePairValues(const string& vec1, const string& vec2) {
    size_t noscansize = primInEdges.size();
    // Full scan
    if (vec1.length() > noscansize) {
        // Flipflop get value from vec
        for(size_t i = 0; i < seqGates.size(); i++)
        {
            Vertex* s = seqGates[i];
            for (size_t j = 0; j < s->outEdges.size(); j++)
            {
                Edge* edg2 = s->outEdges[j];
                string port = s->outPorts[j];
                if (port.compare(".Q") == 0) {
                    edg2->oldVal = (vec1[noscansize + i] == '0' ? 0 : 1);
                    edg2->newVal = (vec2[noscansize + i] == '0' ? 0 : 1);
                }
                else if (port.compare(".QN") == 0) {
                    edg2->oldVal = (vec1[noscansize + i] == '0' ? 1 : 0);
                    edg2->newVal = (vec2[noscansize + i] == '0' ? 1 : 0);
                }
            }
        }
    }
    // no scan
    else {
        //
        for(size_t i = 0; i < seqGates.size(); i++)
        {
            Vertex* s = seqGates[i];
            Edge* edg1 = s->findDataEdge();
            for (size_t j = 0; j < s->outEdges.size(); j++)
            {
                Edge* edg2 = s->outEdges[j];
                string port = s->outPorts[j];
                if (port.compare(".Q") == 0) {
                    edg2->oldVal = edg2->newVal;
                    edg2->newVal = edg1->newVal;
                }
                else if (port.compare(".QN") == 0) {
                    edg2->oldVal = edg2->newVal;
                    edg2->newVal = (edg1->newVal == 0 ? 1 : 0);
                }
            }
        }
    }

    for (size_t i = 0; i < noscansize; i++) {
        primInEdges[i]->oldVal = (vec1[i] == '0') ? 0 : 1;
        primInEdges[i]->newVal = (vec2[i] == '0') ? 0 : 1;
    }

    for (Vertex* v : topoArr)
        v->calculateOutTransition();

    return countTransitions();
}

/*
 * Use vec to simulate the graph
 * both oldVal and newVal will be changed
 * oldVal from newVal
*/ 
int Graph::populateValues(const string& vec) {
    size_t noscansize = primInEdges.size();
    // Full scan
    if (vec.length() > noscansize) {
        // Flipflop get value from vec
        for(size_t i = 0; i < seqGates.size(); i++)
        {
            Vertex* s = seqGates[i];
            for (size_t j = 0; j < s->outEdges.size(); j++)
            {
                Edge* edg2 = s->outEdges[j];
                string port = s->outPorts[j];
                edg2->oldVal = edg2->newVal;
                if (port.compare(".Q") == 0) {
                    edg2->newVal = (vec[noscansize + i] == '0' ? 0 : 1);
                }
                else if (port.compare(".QN") == 0) {
                    edg2->newVal = (vec[noscansize + i] == '0' ? 1 : 0);
                }
            }
        }
    }
    // no scan
    else {
        //
        for(size_t i = 0; i < seqGates.size(); i++)
        {
            Vertex* s = seqGates[i];
            Edge* edg1 = s->findDataEdge();
            for (size_t j = 0; j < s->outEdges.size(); j++)
            {
                Edge* edg2 = s->outEdges[j];
                string port = s->outPorts[j];
                edg2->oldVal = edg2->newVal;
                if (port.compare(".Q") == 0) {
                    edg2->newVal = edg1->newVal;
                }
                else if (port.compare(".QN") == 0) {
                    edg2->newVal = (edg1->newVal == 0 ? 1 : 0);
                }
            }
        }
    }

    for (size_t i = 0; i < noscansize; i++) {
        primInEdges[i]->oldVal = primInEdges[i]->newVal;
        primInEdges[i]->newVal = (vec[i] == '0') ? 0 : 1;
    }

    for (Vertex* v : topoArr)
        v->calculateOutTransition();
    return countTransitions();
}


Edge* Vertex::findDataEdge() {
    for (size_t i = 0; i < inEdges.size(); i++) {
        if (strstr(inPorts[i].c_str(), ".DIN"))
            return inEdges[i];
    }
    return nullptr;
}

Graph* createGraphFromVerilog(string filename, string headername, int& primInCnt)
{
    Graph *g = new Graph(headername);
    // Module file read
    ifstream modulefile(filename);
    vector<string> inVarArr;
    vector<string> outVarArr;
    vector<string> tokens;

    string line;
    int assigncnt = 0;
    if (modulefile.is_open()) {
        // Read module file line by line, delimiter ";"
        while (getline(modulefile, line, ';')) {
            size_t pos = line.find_first_not_of(" \n");
            // Skip blank lines
            if (pos == string::npos) continue;
            line = line.substr(pos);
            // Skip wire, endmodule, tri, comment
            if ((line.compare(0, 4, "wire") == 0) ||
                (line.compare(0, 9, "endmodule") == 0) ||
                (line.compare(0, 3, "tri") == 0) ||
                (line.compare(0, 2, "//") == 0))
                continue;
            for (string::iterator it = line.begin(); it != line.end(); ++it) {
                if (*it == ',' || *it == '(' || *it == ')' || *it == '=')
                    *it = ' ';
            }
            tokenize(line, tokens);

            // Input edges
            if (tokens[0].compare(0, 5, "input") == 0) {
                inVarArr = vector<string>(tokens.begin() + 1, tokens.end());
            }
            // Input edges
            else if (tokens[0].compare(0, 6, "output") == 0) {
                outVarArr = vector<string>(tokens.begin() + 1, tokens.end());
            }
            // Assign edges
            else if (tokens[0].compare(0, 6, "assign") == 0) {
                tokens = {"nb1s1", "AU" + to_string(assigncnt++), ".Q", tokens[1], ".DIN", tokens[2]};
                g->createNewVertexForVlog(tokens);
            }
            // Single gate
            else if (!inVarArr.empty() && !outVarArr.empty()){
                g->createNewVertexForVlog(tokens);
            }
        }
        modulefile.close();
    } else {
        cerr << "Error: Module file " << filename << " is not present, Exiting ... \n";
        delete g;
        exit(EXIT_FAILURE);
    }

    /* find the set of primary input and output edges in the graph */
    /* first break the header into words & store in an arr */
    Edge* edg;
    primInCnt = 0;
    for(string& name : inVarArr)
    {
        edg = g->findEdge(name);
        if (!edg)
              cerr << "Warning: Unconnected primary input edge " << name << endl;
        else if(edg->fromNode)
              cerr << "ERROR: primary edge assertion failure\n";
        else
        {
              edg->bPrimIn = true;
              g->primInEdges.push_back(edg);
              primInCnt++;
        }
    }


    for(string& name : outVarArr)
    {
        edg = g->findEdge(name);
        if (!edg)
              cerr << "Warning: Unconnected primary output edge " << name << endl;
        else
        {
            edg->bPrimOut = true;
            g->primOutEdges.push_back(edg);
        }
    }
    return g;
}

void tokenize(string line, vector<string>& tokens) {
    istringstream iss(line);
    tokens = vector<string>((istream_iterator<string>(iss)),
                            istream_iterator<string>());
}

/* create a vertex of given type - written for generating graph from gate-level verilog mapped to LEDA library 
 * tokens[0] is the type, tokens[1] is the name
*/
Vertex* Graph::createNewVertexForVlog(vector<string>& tokens){
    Edge* edg;
    Vertex* vert = new Vertex(tokens[0], tokens[1]);
    vertexArr.push_back(vert);

    vert->isSeq = strstr(tokens[0].c_str(), "dff");
    if (vert->isSeq) {
        seqGates.push_back(vert);
    } else {
        combGates.push_back(vert);
    }

    // first create the output edges
    for(size_t i = 2; i < tokens.size(); i += 2)
    {
        const char* tstr = tokens[i].c_str();
        if (strstr(tstr, ".Q") || strstr (tstr, ".OUT") || strstr(tstr, ".Out") || strstr(tstr, ".out"))
        {
            edg = createEdge(tokens[i + 1]);
            edg->fromNode = vert;
            vert->outEdges.push_back(edg);
            vert->outPorts.push_back(tokens[i]);
        }
    }

    // now create the input edges
    for(size_t i = 2; i < tokens.size(); i += 2)
    {
        const char* tstr = tokens[i].c_str();
        if (strstr(tstr, ".DIN") || strstr (tstr, ".IN") || strstr(tstr, ".CLK") 
		   || strstr(tstr, ".SIN") || strstr(tstr, ".CIN") || strstr(tstr, ".din") 
		   || strstr(tstr, "clk") || strstr(tstr, ".in") || strstr(tstr, ".AIN")
           || strstr(tstr, ".CLR") || strstr(tstr, ".SET") || strstr(tstr, ".EB") || strstr(tstr, ".BIN"))
        {
            edg = createEdge(tokens[i + 1]);
            edg->toNodes.push_back(vert);
            vert->inEdges.push_back(edg);
            vert->inPorts.push_back(tokens[i]);
        }
    }

    return vert;
}

// Create a new edge with name if not found in g
Edge* Graph::createEdge(string& name){
    Edge* edg = findEdge(name);

    if (edg)
        return edg;

    edg = new Edge(name);
    if (name.compare("1'b0") == 0) {
        edg->name = "GND";
        edg->oldVal = 0;
        edg->newVal = 0;
    }
    else if (name.compare("1'b1") == 0){
        edg->name = "VDD";
        edg->oldVal = 1;
        edg->newVal = 1;
    }
    
    edgArr.push_back(edg);
    return edg;    
}

// Find Edge by name
Edge* Graph::findEdge(const string& name){
    for (size_t i = 0; i < edgArr.size(); i++) {
        if (edgArr[i]->name.compare(name) == 0)
            return edgArr[i];
    }
    return nullptr;
}

// Sort vertices of g in topo order
bool Graph::getTopoOrder(){
    // Assign unique idx to every vertex
    for (size_t i = 0; i < combGates.size(); i++) {
        combGates[i]->idx = i;
    }
    const size_t nNodes = combGates.size();
    for (size_t i = 0; i < seqGates.size(); i++) {
        seqGates[i]->idx = i + nNodes;
    }

    /* topological sort of the vertices */
    vector<int> d(nNodes);
    vector<bool> visited(nNodes, false);
    for (Vertex* v : combGates) {
        int inpCnt = 0;
        for (Edge* edg : v->inEdges) {
            // If input edge is from VDD or GND, skip
            if (!edg->name.compare("VDD") || !edg->name.compare("GND"))
                continue;
            // If input edge is not primary input or seq
            if (!edg->bPrimIn && !(edg->fromNode && strstr(edg->fromNode->type.c_str(), "dff")))
                inpCnt++;
        }
        d[v->idx] = inpCnt;
    }

    Vertex* v = nullptr;
    for (size_t i = 0; i < nNodes; ++i) {
        bool flag = false;
        for (size_t j = 0; j < nNodes; ++j) {
            if (d[j] == 0 && visited[j] == false) {
                visited[j] = true;
                v = combGates[j];
                flag = true;
                break;
            }
        }

        if (!flag) continue;
        topoArr.push_back(v);
        for (Edge* edg : v->outEdges) {
            for (Vertex* v2 : edg->toNodes) {
                if (strstr(v2->type.c_str(), "dff")) continue;
                --d[v2->idx];
            }
        }
    }

    // check topology errors
    for (Vertex* v: combGates) {
        if (d[v->idx]) {
            cerr << "\nERROR: Violation of topological order for vertex " << v->name
                 <<". Number of untraversed input edges: " << d[v->idx] << endl;
            exit(1);
        }
    }

    // Assign level to each vertex
    for (Vertex *v : topoArr) {
        v->level = 0;
        for (Edge *e : v->inEdges) {
            Vertex *t = e->fromNode;
            if (t && (t->level + 1 >= v->level))
                v->level = t->level + 1;
        }
    }
    return true;
}

int Graph::countTransitions() {
    int vectorActivityCnt = 0;
    for (Vertex* s : topoArr) {
        for (Edge *q : s->inEdges) {
            vectorActivityCnt += (q->oldVal != q->newVal) ? 1 : 0;
        }
    }
    for (Edge *q : primOutEdges)
        vectorActivityCnt += (q->oldVal != q->newVal) ? 1 : 0;

    return vectorActivityCnt;
}

// Calculate both old values and new values
int Vertex::calculateOutTransition()
{
    vector<int> oldInVal, newInVal;
    int oVal1, oVal2, nVal1, nVal2;

    /* let's assume that there will be max 2 outputs per cell */
    /* if more, we can pass more oVal to the Vertex::getValue function */

    for(Edge* edg : inEdges)
    {
        oldInVal.push_back(edg->oldVal);
        newInVal.push_back(edg->newVal);
    }

    /* compute the value based on gate type */
    Vertex::getValue(type, oldInVal, oVal1, oVal2);
    Vertex::getValue(type, newInVal, nVal1, nVal2);

    /* set the output values on output edges */
    outEdges[0]->oldVal = oVal1;
    outEdges[0]->newVal = nVal1;
    
    if (outEdges.size() == 2)
    {
        outEdges[1]->oldVal = oVal2;
        outEdges[1]->newVal = nVal2;
    }
    else if (outEdges.size() > 2)
        cerr << "ERROR: out edges > 2 in CalculateOutTransition.\n";
    
    if (oVal1 != nVal1)
        return 1;
    else
        return 0;
}

// Calculate both old values and new values
int Vertex::calculateOutTransitionNewOnly()
{
    vector<int> newInVal;
    int nVal1, nVal2;

    /* let's assume that there will be max 2 outputs per cell */
    /* if more, we can pass more oVal to the Vertex::getValue function */

    for(Edge* edg : inEdges)
    {
        newInVal.push_back(edg->newVal);
    }

    /* compute the value based on gate type */
    Vertex::getValue(type, newInVal, nVal1, nVal2);

    /* set the output values on output edges */
    outEdges[0]->newVal = nVal1;
    
    if (outEdges.size() == 2)
    {
        outEdges[1]->newVal = nVal2;
    }
    else if (outEdges.size() > 2)
        cerr << "ERROR: out edges > 2 in CalculateOutTransition.\n";

    if (outEdges[0]->oldVal != nVal1)
        return 1;
    else
        return 0;
}


void Vertex::getValue(string& type, vector<int>& DIN, int& val1, int& val2)
{
    val1 = -1;
    val2 = -1;

    if(type.compare(0, 3, "nor") == 0)
    {
        val1 = 1;
        for (int i : DIN)
            if (i) val1 = 0;
    }
    else if(type.compare(0, 2, "or") == 0)
    {
        val1 = 0;
        for (int i : DIN)
            if (i) val1 = 1;
    }
    else if(type.compare(0, 3, "nnd") == 0)
    {
        val1 = 0;
        for (int i : DIN)
            if (i == 0) val1 = 1;
    }
    else if(type.compare(0, 3, "and") == 0)
    {
        val1 = 1;
        for (int i : DIN)
            if (i == 0) val1 = 0;
    }
    else if(type.compare(0, 3, "xor") == 0)
    {
        val1 = 1;
        int cnt = 0;
        for (int i : DIN)
            if (i == 1) cnt++;
        if (cnt % 2 == 0)
            val1 = 0;
    }
    else if(type.compare(0, 3, "xnr") == 0)
    {
        val1 = 1;
        int cnt = 0;
        for (int i : DIN)
            if (i == 1) cnt++;
        if (cnt % 2 == 1)
            val1 = 0;
    }
    else if(!type.compare(0, 2, "hi") || !type.compare(0, 2, "i1") || !type.compare(0, 3, "ib1"))
    {
        if(DIN[0] == 0)
            val1 = 1;
        else
            val1 = 0;
    }
    else if(!type.compare(0, 3, "nb1"))
    {
        if(DIN[0] == 0)
            val1 = 0;
        else
            val1 = 1;
    }
    else
        cerr << "ERROR: unknown type of cell: " << type << " in Vertex::getValue.\n";
}

void Graph::populateAfterFlipping(Edge* edg)
{
    Vertex *from = edg->fromNode;
    auto it = topoArr.begin();
    while (it != topoArr.end() && (*it) != from) ++it;
    if (it == topoArr.end()) {
        cerr << "Error in finding from node to " << edg->name << endl;
        return;
    }
    edg->newVal = 1 - edg->newVal;
    while ((++it) != topoArr.end()) {
        (*it)->calculateOutTransitionNewOnly();
    }
}
