/*
 * Authors: Yangdi Lyu
 * Jan 2018
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stack>
#include <unordered_set>
#include "graph.h"
#include "traverse.h"
#include "genetic.h"

config_t config;
static Edge* Payload;

// Simulate with filename
// Get info about the rare nodes
void generateStat(Graph *g, const string& filename, float TRIGTHRESH) {
    ifstream fvec(filename);
    if (!fvec.is_open()) {
        cerr << "File for generate stat not present: " << filename << ". Exiting ...\n";
        exit(EXIT_FAILURE);
    }
    int nPatterns = 0,  activityCnt = 0, zerocnt = 0, onecnt = 0;

    cout << "\nAnalyzing internal node activities for random vectors from " << filename << " ...\n";

    // Initialize all signals to zero
    // TODO: Use primary input + out edge of vertex to avoid assign to edges not connected, e.g., VDD
    for (Edge *q : g->edgArr)
    {
        q->oldVal = q->newVal = q->signal0Prob = q->signal1Prob = 0;
        q->activity01 = q->activity10 = q->lowProbVal = 0;
    }
    for(Vertex* s: g->topoArr)
    {
        s->nInput00 = s->nInput01 = s->nInput10 = s->nInput11 = 0;
        if (s->outEdges.size()) {
            for (Edge* q : s->outEdges) {
                q->oldVal = q->newVal = q->signal0Prob = q->signal1Prob = 0;
                q->activity01 = q->activity10 = q->lowProbVal = 0;
            }
        }
    }

    // Read test vectors from file
    string vec1;
    while(fvec.good())
    {
        fvec >> vec1;
        if (!fvec.good())
            break;

        // apply the test pattern to PI's
        g->populateValues(vec1);
        activityCnt += g->countTransitions();
        for(Edge* q : g->primInEdges)
        {
            q->signal0Prob = (q->newVal) ? q->signal0Prob : q->signal0Prob + 1;
            q->signal1Prob = (q->newVal) ? q->signal1Prob + 1 : q->signal1Prob;
            q->activity01 = (q->oldVal == 0 && q->newVal == 1) ? q->activity01 + 1 : q->activity01;
            q->activity10 = (q->oldVal == 1 && q->newVal == 0) ? q->activity10 + 1 : q->activity10;
        }

        for(Vertex* s : g->topoArr)
        {
            Edge* edg1 = s->inEdges[0];
            Edge* edg2 = (s->inEdges.size() > 1) ? s->inEdges[1] : nullptr;
            s->nInput00 = (!edg1->oldVal && ((edg2 && !edg2->oldVal) || !edg2)) ? s->nInput00 + 1 : s->nInput00;
            s->nInput01 = (!edg1->oldVal && edg2 && edg2->oldVal) ? s->nInput01 + 1 : s->nInput01;
            s->nInput10 = (edg1->oldVal && edg2 && !edg2->oldVal) ? s->nInput10 + 1 : s->nInput10;
            s->nInput11 = (edg1->oldVal && ((edg2 && edg2->oldVal) || !edg2)) ? s->nInput11 + 1 : s->nInput11;
            
            s->nInput00 = (!edg1->newVal && ((edg2 && !edg2->newVal) || !edg2)) ? s->nInput00 + 1 : s->nInput00;
            s->nInput01 = (!edg1->newVal && edg2 && (edg2 && edg2->newVal)) ? s->nInput01 + 1 : s->nInput01;
            s->nInput10 = (edg1->newVal && edg2 && !edg2->newVal) ? s->nInput10 + 1 : s->nInput10;
            s->nInput11 = (edg1->newVal && ((edg2 && edg2->newVal) || !edg2)) ? s->nInput11 + 1 : s->nInput11;
        }

        for(Vertex* s : g->topoArr)
        {
            if (s->outEdges.size())
            {
                Edge *q = s->outEdges[0];
                q->signal0Prob = (q->newVal) ? q->signal0Prob : q->signal0Prob + 1;
                q->signal1Prob = (q->newVal) ? q->signal1Prob + 1 : q->signal1Prob;
                q->activity01 = (q->oldVal == 0 && q->newVal == 1) ? q->activity01 + 1 : q->activity01;
                q->activity10 = (q->oldVal == 1 && q->newVal == 0) ? q->activity10 + 1 : q->activity10;
            }
        }
		  
        nPatterns++;
    }

    fvec.close();

    // Simulation done, computing the probability of each value
    for(Edge* q : g->edgArr)
    {
        q->signal0Prob = (float)q->signal0Prob/nPatterns;
        q->signal1Prob = (float)q->signal1Prob/nPatterns;
        q->activity01 = (float)q->activity01/nPatterns;
        q->activity10 = (float)q->activity10/nPatterns;
    }

    //  -------------------------------- rare value is 0 --------------------------------------//
    for(Edge* q : g->primInEdges)
        if (q->signal0Prob <= TRIGTHRESH) {
            g->lowprobEdges.push_back(q);
            q->lowProbVal = 0;
        }

    for(Vertex* s : g->topoArr)
    {
        if (s->outEdges.size())
        {
            Edge* q = s->outEdges[0];
            if (q->signal0Prob <= TRIGTHRESH && q->name.compare("VDD"))
            {
                zerocnt++;
                g->lowprobEdges.push_back(q);
                q->lowProbVal = 0;
            }
        }
    }

    //  -------------------------------- rare value is 1 --------------------------------------//
    for(Edge* q : g->primInEdges)
        if (q->signal1Prob <= TRIGTHRESH) {
            g->lowprobEdges.push_back(q);
            q->lowProbVal = 1;
        }

    for(Vertex* s : g->topoArr)
    {
        if (s->outEdges.size())
        {
            Edge* q = s->outEdges[0];
            if (q->signal1Prob <= TRIGTHRESH && q->name.compare("GND"))
            {
                onecnt++;
                g->lowprobEdges.push_back(q);
                q->lowProbVal = 1;
            }
        }
    }
    
	cout << "\n\nNumber of patterns: " << nPatterns << endl;
    cout << "\n\nGraph Statistics ...\n";
    cout << "Graph vertices: " << g->vertexArr.size() << endl;
    cout << "Number of flip-flops: " << g->seqGates.size() << endl;
    cout << "Number of Edges with Signal Probability below Threshold (" << TRIGTHRESH << "): " << onecnt + zerocnt <<  endl;
    cout << "Graph Edges: " << g->edgArr.size() << endl;
    cout << "Average activity per node: " << (float)activityCnt/(float)(nPatterns * g->topoArr.size())
         << ". ActivityCnt per test: " << (float)activityCnt/(float)nPatterns << endl;;
         
    cout << "\nINFO: Number of Edges with Signal Probability below Threshold (" << TRIGTHRESH << "): " << onecnt + zerocnt <<  endl;

}

// Generating N-detect test vectors
// Compute coverage of random test vectors and N-detect test vectors
void Ndetect(Graph *G, const string& randFileName, const string& firstvecfile)
{
    if (access(firstvecfile.c_str(), F_OK ) != -1 ) {
        cout << "N-detect file: " << firstvecfile << " already exist. Skip generating...\n";
        return;
    }
    // Use multimap to sort random tests
    multimap<int, string> testset;
    sortRandomTestset(G, randFileName, testset);
    // Here is MERO
    generateGoodVectors(G, testset, firstvecfile);
}

// Simulate test vectors from filename
// Sort based on the number of hits in rare nodes
void sortRandomTestset(Graph* G, const string& filename, multimap<int, string>& testset) {
    ifstream fvec(filename);
    if (!fvec.is_open()) {
        cerr << "Vector file not present: " << filename << ". Exiting ...\n";
        exit(EXIT_FAILURE);
    }
    // Initialize the graph
    for (Edge* q : G->primInEdges)
        q->oldVal = q->newVal = 0;
    for (Vertex* v : G->vertexArr)
        for (Edge*q : v->outEdges)
            q->oldVal = q->newVal = 0;

    // Read vector file line by line
    string line;
    
    while (fvec.good()) {
        fvec >> line;
        if (!fvec.good())	break;

        int nhits = 0;
        G->populateValues(line);
        for (Edge* edg : G->lowprobEdges)
            nhits += (edg->newVal == edg->lowProbVal) ? 1 : 0;

        if (nhits > 0)
            testset.insert(make_pair(nhits, line));
    }
    fvec.close();
}

// Simulate with test vectors from filename
// triggers and payload already read from file
void triggerSim(Graph* G, vector<vector<Edge*>>& ntriggers, vector<Edge*>& npayloads, const string& filename) {
    int trigcoveredCnt = 0;
    int trojcoveredCnt = 0;
    ifstream fvec(filename);
    if (!fvec.is_open()) {
        cerr << "Vector file not present: " << filename << ". Exiting ...\n";
        exit(EXIT_FAILURE);
    }
    // Initialize the graph
    for (Edge* q : G->primInEdges)
        q->oldVal = q->newVal = 0;
    for (Vertex* v : G->vertexArr)
        for (Edge*q : v->outEdges)
            q->oldVal = q->newVal = 0;

    // Read vector file line by line
    string line;
    
    vector<int> trigcoveredList(npayloads.size(), 0);
    vector<bool> trojcoveredList(npayloads.size(), false);
    
    while (fvec.good()) {
        fvec >> line;
        if (!fvec.good())	break;

        int nhits = 0;
        G->populateValues(line);
        for (Edge* edg : G->lowprobEdges)
            nhits += (edg->newVal == edg->lowProbVal) ? 1 : 0;

        // Check for trigger conditions
        for (size_t i = 0; i < npayloads.size(); ++i){
            if (trojcoveredList[i])	continue;
            bool covered = true;
            for (Edge* edg : ntriggers[i]) {
                if (edg->newVal != edg->lowProbVal){
                    covered = false;
                    break;
                }
            }
            if (covered) {
                if (trigcoveredList[i] < config.SEQSTATES) {
                    ++trigcoveredList[i];
                    if (trigcoveredList[i] == config.SEQSTATES) {
                        trigcoveredCnt++;
                    }
                }
                
                // Without checking if the payload propogate to primary output
                if (trigcoveredList[i] == config.SEQSTATES) {
                    trojcoveredList[i] = true;
                    trojcoveredCnt++;
                }
            }
        }
    }
    fvec.close();

    cout << "\nResult (" << filename << ") with Random Vectors: Suspects: " << G->lowprobEdges.size() << " (Out of: "
    << G->edgArr.size() <<  ") Instances: " << ntriggers.size() << " Triggers Covered: " << trigcoveredCnt <<" Trojans Covered: "
    << trojcoveredCnt << " Trigger Cov: " << 100 * (float)trigcoveredCnt / ntriggers.size() << " Trojan Cov:"
    << 100 * ((float)trojcoveredCnt/ntriggers.size()) << endl;
}

// N-detect test generation
void generateGoodVectors(Graph* G, multimap<int, string>& testset, const string& firstvecfile){
	ofstream fvec(firstvecfile);

    cout << "INFO: Generating N-Detect vectors ...\n";
    string vec(testset.begin()->second.length(), '0');
    G->populateValues(vec);

    vector<int> trigList(G->lowprobEdges.size(), 0);

    for (auto it = testset.rbegin(); it != testset.rend(); it++) {
        // Flip bits of test vectors to increase coverage
        string vec2 = modifyVectorToImproveCoverage(G, trigList, it->second);
        G->populateOnlyNewValues(vec2);

        int nhits = countHitsWithTriglist(G, trigList, true);
        if (nhits)
            fvec << vec2 << endl;
        
        // Let the new values in old values
        G->populateValues(vec2);
    }
    fvec.close();
}

// Modify test vector vec to increase coverage
string modifyVectorToImproveCoverage(Graph* G, vector<int>& trigList, const string& vec){
    G->populateOnlyNewValues(vec);
    // Count the number of hits without updating
    int besthit = countHitsWithTriglist(G, trigList, false);

    // ROUNDS of one-bit flipping
    int ROUNDS = 2;
    string bestvec = vec;
    for (int r = 0; r < ROUNDS; r++) {
        for (size_t j = 0; j < bestvec.length(); j++) {
            bestvec[j] = (bestvec[j] == '0') ? '1' : '0';
            G->populateOnlyNewValues(bestvec);

            int nhits = countHitsWithTriglist(G, trigList, false);
            if (nhits <= besthit)
                bestvec[j] = (bestvec[j] == '0') ? '1' : '0';
            else
                besthit = nhits;
        }
    }

    // One round of two-bit flipping
    for (size_t j = 0; j + 1 < bestvec.length(); j++) {
        bestvec[j] = (bestvec[j] == '0') ? '1' : '0';
        bestvec[j + 1] = (bestvec[j + 1] == '0') ? '1' : '0';
        G->populateOnlyNewValues(bestvec);

        int nhits = countHitsWithTriglist(G, trigList, false);
        if (nhits <= besthit){
            bestvec[j] = (bestvec[j] == '0') ? '1' : '0';
            bestvec[j + 1] = (bestvec[j + 1] == '0') ? '1' : '0';
        }
        else
            besthit = nhits;
    }

    return bestvec;

}

// Count the number of rare nodes hit
int countHitsWithTriglist(Graph* G, vector<int>& trigList, bool update) {
    int nhits = 0;
    // If update, inc = 1
    int inc = update ? 1 : 0;
    for (size_t i = 0; i < G->lowprobEdges.size(); i++) {
        Edge* q = G->lowprobEdges[i];
        if (trigList[i] < config.NDETECT && (q->newVal == q->lowProbVal)) {
            nhits++;
            trigList[i] += inc;
        }
    }
    return nhits;
}

// Read trojans from file
// Trojans are defined as [Trigger1, Trigger2, ..., TriggerN, Payload] in trojanfile 
void readTrojansFromFile(Graph* G, const string& trojanfilename, vector<vector<Edge*>>& ntriggers, vector<Edge*>& npayloads){
    ifstream ftroj(trojanfilename);
    if (!ftroj.is_open()) {
        cerr << "Trojan file not present: " << trojanfilename << ". Exiting ...\n";
        exit(EXIT_FAILURE);
    }
    cout << "Reading Trojan from file: " << trojanfilename << endl;
    string line;
    while (getline(ftroj, line)) {
        istringstream iss(line);
        vector<Edge*> triggers;
        string edgname;   
        while (iss.good()) {
            iss >> edgname;
            triggers.push_back(G->findEdge(edgname));
        }
        npayloads.push_back(triggers.back());
        triggers.pop_back();
        ntriggers.push_back(triggers);
    }
    ftroj.close();
    Payload = npayloads[0];
}

// Run GA to generate pairs for each vector from vecfile
void genSuccVectorsGA(Graph* G, const string& firstvecfile, const string& secondvecfile, vector<Edge*>& npayloads)
{
    string trojanfilename = string("./Trojans_theta_") + config.TRIGTHRESH + "/" + config.moduleName +  ".v_" + to_string(config.NTRIG) + ".trojans_" + to_string(config.TRIGGER_INSTANCES);

    if (access(secondvecfile.c_str(), F_OK) != -1) {
        cout << "GA file: " << secondvecfile << " already exist. Skip generating...\n";
        return;
    }
	
    // Input file
    ifstream fvec(firstvecfile);
    ofstream fout(secondvecfile);
    if (!fvec.is_open()) {
        cerr << "File not present: " << firstvecfile << ". Exiting ...\n";
        exit(EXIT_FAILURE);
    }
    //  Read all vecs
    string vec, pair;
    while(fvec.good())
    {
        fvec >> vec;
        if(!fvec.good())  break; 
        // Run GA to generate pair for vec
        geneOneVectorGA(G, vec, pair);
        fout << pair << endl;
    }
    fvec.close();
    fout.close(); 
}

// Run BFS to generate pairs for each vector from vecfile
void genSuccVectorsBFS(Graph* G, const string& firstvecfile, const string& secondvecfile)
{
    if ( access(secondvecfile.c_str(), F_OK ) != -1 ) {
        cout << "BFS file: " << secondvecfile << " already exist. Skip generating...\n";
        return;
    }
	
    // Input file
    ifstream fvec(firstvecfile);
    ofstream fout(secondvecfile);
    if (!fvec.is_open()) {
        cerr << "File not present: " << firstvecfile << ". Exiting ...\n";
        exit(EXIT_FAILURE);
    }
    //  Read all vecs
    string vec, pair;
    while(fvec.good())
    {
        fvec >> vec;
        if(!fvec.good())  break; 

        fout << geneOneVectorBFS(G, vec) << endl;
    }
    fvec.close();
    fout.close(); 
}

// Generate one pair for vec1
// Here we use two levels of traverse
// TODO: user define the number of levels
string geneOneVectorBFS(Graph* G, const string& vec1){
    string vec2 = vec1, bestvec = vec1;
    double bestfit = -1.0, singlefit = -1.0;
    // First level
    for (size_t i = 0; i < vec1.length(); i++) {
        vec2[i] = (vec2[i] == '0') ? '1' : '0';
        singlefit = switchSimPatternOnePair(G, vec1, vec2);
        if (singlefit > bestfit) {
            bestfit = singlefit;
            bestvec = vec2;
        }
        vec2[i] = (vec2[i] == '0') ? '1' : '0';
    }

    // Second level
    for (size_t i = 0; i < vec1.length(); i++) {
        for (size_t j = i + 1; j < vec1.length(); j++) {
            vec2[i] = (vec2[i] == '0') ? '1' : '0';
            vec2[j] = (vec2[j] == '0') ? '1' : '0';
            singlefit = switchSimPatternOnePair(G, vec1, vec2);
            if (singlefit > bestfit) {
                bestfit = singlefit;
                bestvec = vec2;
            }
            vec2[i] = (vec2[i] == '0') ? '1' : '0';
            vec2[j] = (vec2[j] == '0') ? '1' : '0';
        }
    }
    return bestvec;
}

// Read from trojan files, create a module with a trojan for each
void createTrojans(vector<Graph*>& Gt, Graph* G, int instances)
{
    string headername = G->header;
    char trojanfilename[100];
    sprintf(trojanfilename, "./Trojans_theta_%s/%s.v_%d.trojans_%d", config.TRIGTHRESH.c_str(), G->header.c_str(), config.NTRIG, instances);
    ifstream ftrj(trojanfilename);

    if (!ftrj.is_open()) {
        cerr << "Trojan file not present: " << trojanfilename << " Exiting ...\n";
        exit(EXIT_FAILURE);
    }

    cout << "Creating trojans with file " << trojanfilename << endl;

    vector<string> triggerIns(config.NTRIG);
    string line, payload;
    while (ftrj.good()) {
        for (int i = 0; i < config.NTRIG; i++) {
            ftrj >> line;
            triggerIns[i] = line;
        }
        ftrj >> payload;
        if (!ftrj.good())   break;

        // Create a new graph
        int primInCnt;
        Graph* gh = createGraphFromVerilog(config.inputFile, headername, primInCnt);
        for (Edge* edg : G->lowprobEdges)
        {
            gh->lowprobEdges.push_back(gh->findEdge(edg->name));
            Edge* tedg = gh->findEdge(edg->name); 
            tedg->lowProbVal = edg->lowProbVal; 
            gh->lowprobEdges.push_back(tedg);
        }

        // Insert trojan
        insertTrojanNew(gh, triggerIns, payload);
        Gt.push_back(gh);
        gh->getTopoOrder();
    }
    
    ftrj.close();
}

// read in PI vectors and simulate the G and Gt(W/ trojan) side by side
//       get the switching difference and relative switching 
// vec1: holds the first vector of the pair;  vec2: holds the second vector of the pair;
double switchSimPatternOnePair(Graph* G, const string& vec1, const string& vec2)
{
    int sw1 = G->populatePairValues(vec1, vec2);
    if ((vec1.compare(vec2) == 0) || (sw1 == 0))	return 0;
    // simple fitness function for DATE paper
    // return (countRareSwitch(G)) / (double)sw1;
    G->populateAfterFlipping(Payload);
    int sw2 = G->countTransitions();
    int switchcnt = (sw1 > sw2) ? (sw1 - sw2) : (sw2 - sw1);
    return  countRareSwitch(G) * ((switchcnt > 0.3 * sw1) ? 0.3 : (switchcnt / (double)sw1));
}

// Count Rare Switch
int countRareSwitch(Graph* G){
    int count = 0;
    for (Edge* e : G->lowprobEdges) {
        if (e->oldVal != e->newVal)
            count++;
    }
    return count;
}

// Read in PI vectors and simulate the G and Gt(W/ trojan) side by side
// get the switching difference and relative switching 
// prevfile: holds the first vector of the pair;  vectfile: holds the second vector of the pair;
void switchSimPatternPair(Graph* G, const string& firstvecfile, const string& secondvecfile)
{
    vector<Graph*> Gt;
    createTrojans(Gt, G, config.TRIGGER_INSTANCES);
    cout << "Created " << config.TRIGGER_INSTANCES << " Trojans\n";

    ifstream fvec1(firstvecfile);
    ifstream fvec2(secondvecfile);
    if (!fvec1.is_open() || !fvec2.is_open()) {
        cerr << "Patterns not present: " << firstvecfile << " , or " << secondvecfile << ". Exiting...\n";
        exit(EXIT_FAILURE);
    }
    cout << "Start simulating with pairs from " << firstvecfile << " (1st vec) and " << secondvecfile << " (2nd vec).\n";

    int total_orig_switch = 0;
    vector<int> eff_orig_switch(config.TRIGGER_INSTANCES, -1);
    vector<double> rel_switch(config.TRIGGER_INSTANCES, -1);
    string vec1, vec2;
    int nPatterns = 0;
    while (fvec1.good() && fvec2.good()) {
        fvec1 >> vec1;
        fvec2 >> vec2;
        if(!fvec1.good() || !fvec2.good())  break;

        int sw1 = G->populatePairValues(vec1, vec2);
        total_orig_switch += sw1;
	    if ((vec1.compare(vec2) == 0) || (sw1 == 0))	continue;

        for (int i = 0; i < config.TRIGGER_INSTANCES; ++i) {
            int sw2 = Gt[i]->populatePairValues(vec1, vec2);
            // For each Trojan, keep the most anomalous behavior
            double rs = (double)abs(sw1-sw2)/ (double)sw1;
            if (rs > rel_switch[i]) {
                rel_switch[i] = rs;
                eff_orig_switch[i] = sw1;
            }
        }
        ++nPatterns;
    }
    fvec1.close();
    fvec2.close();
    // Average relative switch over all trojans
    if (config.debug) {
        ofstream fdebug(string("Debug/") + config.moduleName + ".debug");
        fdebug << "Effective orig switch \t relative switch \n";
        for (int i = 0; i < config.TRIGGER_INSTANCES; ++i)
            fdebug << eff_orig_switch[i] << "\t" << rel_switch[i] << endl;
        fdebug.close();
    }
    double avg_rel_switch = std::accumulate(rel_switch.begin(), rel_switch.end(), 0.0) / rel_switch.size();
    cout << "\nResult: " << firstvecfile << "->" << secondvecfile
         << ", orig_switch, rel_switch\n" <<  "\nResult: " << firstvecfile << "->" << secondvecfile
         << ", " << (total_orig_switch / (double) nPatterns) << ", "<< avg_rel_switch << "\n";  
}

// Insert a new trojan into graph G
void insertTrojanNew(Graph* G, vector<string>& triggerIns, const string& payload)
{
    Edge* edg = G->findEdge(payload);
    // Check if payload exists
    if (edg == nullptr) {
        cerr << "Payload " << payload << " is not found. Exiting ... \n";
        exit(EXIT_FAILURE);
    }
    string oldname = edg->name, newname = edg->name + "_ORG";
    Vertex* fromv = edg->fromNode;
    if (fromv == nullptr) {
        cerr << "Payload " << payload << " cannot be primary input. Exiting ... \n";
        exit(EXIT_FAILURE);
    }

    // Create new edge as the out edge for fromv
    Edge* nedg = G->createEdge(newname);
    fromv->outEdges.clear();
    fromv->outEdges.push_back(nedg);
    nedg->fromNode = fromv;

    // Create trigger gates
    int trjGateCnt = 0;
    queue<string> current;
    vector<string> tokens;
    for (vector<string>::iterator it = triggerIns.begin(); it != triggerIns.end(); ++it) {
        // Check if trigger exists
        Edge *q = G->findEdge(*it);
        if (q == nullptr) {
            cerr << "Trigger net " << *it << " is not found. Exiting ... \n";
            exit(EXIT_FAILURE);
        }
        // If payload is also trigger
        if (payload.compare(*it) == 0)
            *it = newname;

        if (q->lowProbVal == 1) {
            current.push(*it);
        } else {
            // Insert an inverter
            tokens = {"hi1s1", "Trj_" + to_string(trjGateCnt++), ".Q", *it + "_inv", ".DIN", *it};
            G->createNewVertexForVlog(tokens);
            current.push(*it + "_inv");
        }
    }

    // and2s1 gates
    tokens = {"and2s1", "", ".Q", "", ".DIN1", "", ".DIN2", ""};
    while (current.size() > 1) {
        tokens[1] = "Trj_" + to_string(trjGateCnt);
        tokens[3] = "T" + to_string(trjGateCnt++) + "_1";
        tokens[5] = current.front();
        current.pop();
        tokens[7] = current.front();
        current.pop();
        G->createNewVertexForVlog(tokens);
        current.push(tokens[3]);
    }

    // xor2s1 gate (payload gate)
    tokens = {"xor2s1", "Trj_" + to_string(trjGateCnt), ".Q", oldname, ".DIN1", newname, ".DIN2", current.front()};
    G->createNewVertexForVlog(tokens);
}

// Generate random trojans
// Need TMAX to validate each trojan
void genRandTrojans(Graph* G){
    char trojanfilename[100];
    sprintf(trojanfilename, "./Trojans_theta_%s/%s.v_%d.payloads_%d", config.TRIGTHRESH.c_str(), G->header.c_str(), config.NTRIG, config.TROJAN_INSTANCES / 10);

    if ( access(trojanfilename, F_OK ) != -1 ) {
        cout << "Payload file: " << trojanfilename << " already exist. Skip generating...\n";
        return;
    }

    cout << "INFO: Creating random trojan instances ...\n";

    srand(4321);

    // generate random payloads
    char payloadfile[100];
    sprintf(payloadfile, "./Trojans_theta_%s/%s.v_%d.payloads_%d", config.TRIGTHRESH.c_str(), G->header.c_str(), config.NTRIG, config.TROJAN_INSTANCES / 10);

    vector<int> payloads;
    ofstream fpayloads(payloadfile);
    int total = G->topoArr.size() / 2;
    for (int i = 0; i < config.TRIGGER_INSTANCES / 10; i++) {
        int pyld = rand() % (total--) + G->topoArr.size() / 2;
        while(find(payloads.begin(), payloads.end(), pyld) != payloads.end()) {
            pyld++;
        }
        payloads.push_back(pyld);
        fpayloads << G->topoArr[pyld]->name << endl;
    }
    fpayloads.close();

    // Create Trojans
    ofstream ftrj("trojans.txt");
    ftrj << "***                        Trigger edges                    ***    Payload      ***\n";
    ftrj << "---------------------------------------------------------------------------\n";
    for (int i = 0; i < 1000000; i++) {
        ftrj << i << ". ";
        int nlow = G->lowprobEdges.size();
        vector<int> triggers;
        // Randomly select triggers
        for (int j = 0; j < config.NTRIG; ++j) {
            int t = rand() % (nlow--);
            while(find(triggers.begin(), triggers.end(), t) != triggers.end()) {
                t++;
            }
            triggers.push_back(t);
            Edge* e = G->lowprobEdges[t];
            ftrj << e->name << "(val: " << e->lowProbVal << ", idx: " << e->fromNode->idx << ")    ";
        }
        // Randomly select payloads
        Edge* e = G->topoArr[payloads[rand() % payloads.size()]]->outEdges[0];
        ftrj << "    " << e->name << "(idx: " << ((e->fromNode) ? e->fromNode->idx : 0) << ")\n";
    }
    ftrj << "\n-----------------------------------------------------------------------------\n";
    ftrj.close();

    // Run TMAX to validate
    char cmdstr[100];
    sprintf(cmdstr,"./tj_justify_triggers %s.v %s %s X X 0 X", config.moduleName.c_str(), config.moduleName.c_str(), config.cktType.c_str());
    // Execute the TCL code
    cout << cmdstr << endl;
    // system(cmdstr);
}

