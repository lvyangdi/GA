#ifndef __TRAVERSE_H__
#define __TRAVERSE_H__

#include <vector>
#include <string>
#include <map>
#include <unordered_set>
using namespace std;

class Graph;
class Edge;

struct config_t{
    bool fullScan;
    string inputFile;
    int NTRIG;
    int nMaxPatterns;
    double inputActivity;
    int NDETECT;
    int SEQSTATES;
    int NTESTPTS;
    string TRIGTHRESH;
    int TRIGGER_INSTANCES;
    int TROJAN_INSTANCES;
    string moduleName;
    string cktType;
    string isScan;
    bool debug;
};

void generateStat(Graph *g, const string& filename, float TRIGTHRESH);

void sortRandomTestset(Graph* G, const string& randFileName, multimap<int, string>& testset);

void triggerSim(Graph* G, vector<vector<Edge*>>& ntriggers, vector<Edge*>& npayloads, const string& filename);

void Ndetect(Graph *G, const string& randomfile, const string& firstvecfile);

void genSuccVectorsGA(Graph* G, const string& firstvecfile, const string& secondvecfile, vector<Edge*>& npayloads);

void genSuccVectorsBFS(Graph* G, const string& firstvecfile, const string& secondvecfile);

string geneOneVectorBFS(Graph* G, const string& vec1);

void createTrojans(vector<Graph*>& Gt, Graph* G, int instances);

void readTrojansFromFile(Graph* G, const string& trojanfilename, vector<vector<Edge*>>& ntriggers, vector<Edge*>& npayloads);

double switchSimPatternOnePair(Graph* G, const string& vec1, const string& vec2);

void switchSimPatternPair(Graph* G, const string& firstvecfile, const string& secondvecfile);

void insertTrojanNew(Graph* G, vector<string>& triggerIns, const string& payload);

void generateGoodVectors(Graph* G, multimap<int, string>& testset, const string& meroFileName);

string modifyVectorToImproveCoverage(Graph* G, vector<int>& trigList, const string& vec);

int countHitsWithTriglist(Graph* G, vector<int>& trigList, bool update);

int countRareSwitch(Graph* G);

void genRandTrojans(Graph* G);

#endif
