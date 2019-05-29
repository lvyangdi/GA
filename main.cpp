/*
 * Authors: Yangdi Lyu
 * Jan 2018
 */

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <cstdlib>
#include <random>
#include "graph.h"
#include "traverse.h"
using namespace std;

extern config_t config;

void autoconfig(int argc, char* argv[]);
void genRandVectors(string filename, size_t length, unsigned int patterns);
void readLowNodes(Graph *G, std::string nodefilename);
void dumpLowNodes(Graph *G, std::string nodefilename);

int main(int argc, char* argv[]) {
    if(argc < 2)
    {
        cout << "Usage: tcov <Verilog input file> [# of vectors] [trigger threshold] [# of trigger points] [# of trigger/trojan instances] [module name] [comb|seq] [N] [# of seq Trojan states] [# of observable test points] [input activity] [-fullscan | -noscan] [list of trigger nets] [payload net] \n";
        return 1;
    }

    autoconfig(argc, argv);
    const char* testFolder = "./MERO_tests";
    const char* debugFolder = "./Debug";
    if ((mkdir(testFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) && (errno != EEXIST)) {
        std::cerr << "Error creating folder for test vectors" << std::endl;
        return -1;
    }
    if ((mkdir(debugFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) && (errno != EEXIST)) {
        std::cerr << "Error creating folder for test vectors" << std::endl;
        return -1;
    }

    cout << "\nCreating the graph ...\n";
    int nPrimIn;
    Graph *G = createGraphFromVerilog(config.inputFile, config.moduleName, nPrimIn);

    // Generate random test vectors
    string randFileName = string(testFolder) + "/rand" + to_string(config.nMaxPatterns) + "_PI_" + G->header + ".txt";
    string nodeFileName = string(debugFolder) + "/" + config.moduleName + "_nodes_" + config.TRIGTHRESH + ".txt";
    string firstVecFile = string(testFolder) + "/test_" + config.moduleName + "_N" + to_string(config.NDETECT) + "_" + config.TRIGTHRESH + ".txt.first";
    string secondVecFile = string(testFolder) + "/test_" + config.moduleName + "_N" + to_string(config.NDETECT) + "_" + config.TRIGTHRESH + ".txt.second";
    string trojanfilename = string("./Trojans_theta_") + config.TRIGTHRESH + "/" + config.moduleName +  ".v_" + to_string(config.NTRIG) + ".trojans_" + to_string(config.TRIGGER_INSTANCES);
    
    if (config.nMaxPatterns > 0) {
        int length = G->primInEdges.size() + ((config.fullScan) ? G->seqGates.size() : 0);
        genRandVectors(randFileName, length, config.nMaxPatterns);
    }
    
    vector<Edge*> npayloads;
    vector<vector<Edge*>> ntriggers;
    // Read trojans from file
    readTrojansFromFile(G, trojanfilename, ntriggers, npayloads);

    cout << "Getting topologically sorted graph vertices ...\n";
    G->getTopoOrder();
    cout << "\nINFO: Number of vertices: " << G->vertexArr.size()
        << " Number of edges: " << G->edgArr.size() << ".\n";

    // Simulate with random test vectors
    struct stat buffer;
    if (stat(nodeFileName.c_str(), &buffer) == 0) {
        cout << "Low prob nodes file exist. Skip simulation.\n";
        cout << "Reading low prob nodes from file " << nodeFileName << std::endl;
        readLowNodes(G, nodeFileName);
    } else {
        cout << "Start simulation with random vectors.\n";
        generateStat(G, randFileName, atof(config.TRIGTHRESH.c_str()));
        cout << "Dumping low prob nodes to file " << nodeFileName << std::endl;
        dumpLowNodes(G, nodeFileName);
    }

    clock_t start_t, finish_t; 
    
    start_t = clock();
    // Using N-detect to generate preceding vectors
    cout << "\nGenerating first vector in file " << firstVecFile <<" ...\n";
    Ndetect(G, randFileName, firstVecFile);
    finish_t = clock();
    double duration = ((double)(finish_t - start_t)) / CLOCKS_PER_SEC;
    cout << "Result: generate N-detect vectors Runtime: " << duration << " seconds\n";
    
    // Simulate and compute coverage
    triggerSim(G, ntriggers, npayloads, randFileName);
    triggerSim(G, ntriggers, npayloads, firstVecFile);

    start_t = clock();
    // Given preceding vectors, generate succeeding vectors to form pairs
    genSuccVectorsGA(G, firstVecFile, secondVecFile, npayloads);
    finish_t = clock();
    duration = ((double)(finish_t - start_t)) / CLOCKS_PER_SEC;
    cout << "Result: generate succ vectors Runtime: " << duration << " seconds\n";
    switchSimPatternPair(G, firstVecFile, secondVecFile);

    return 0;
}

/*
 * Config the parameters
*/
void autoconfig(int argc, char* argv[]) {
    config.fullScan = false;
    config.NTRIG = 4;
    config.nMaxPatterns = 100000;
    config.inputActivity = 0.50;
    config.NDETECT = 1000;
    config.SEQSTATES = 1;
    config.NTESTPTS = 0;
    config.TRIGTHRESH = "0.1";
    config.TRIGGER_INSTANCES = 1000;
    config.TROJAN_INSTANCES = 1000;
    config.moduleName = "c2670";
    config.cktType = "comb";
    config.inputFile = argv[1];
    config.isScan = "-genePatternPairGA";
    config.debug = true;
    if (argc >= 3)
        config.nMaxPatterns = atoi(argv[2]);
    if (argc >= 4)
        config.TRIGTHRESH = argv[3];
    if (argc >= 5)
        config.NTRIG = atoi(argv[4]);
    if (argc >= 6)
        config.TRIGGER_INSTANCES = atoi(argv[5]);
    if (argc >= 7)
        config.moduleName = ((char **)argv)[6];
    if (argc >= 8)
        config.cktType = ((char **)argv)[7];
    if (argc >= 9)
        config.NDETECT = atoi(argv[8]);
    if (argc >= 10)
        config.SEQSTATES = atoi(argv[9]);
    if (argc >= 11)
        config.NTESTPTS = atoi(argv[10]);
    if (argc >= 12)
        config.inputActivity = atof(argv[11]);
    if (argc >= 13)
        config.isScan = ((char**)argv)[12];

    if (config.cktType.compare("comb"))
        config.fullScan = true;
    config.SEQSTATES = 1;
}

// Generate random test vectors
void genRandVectors(std::string filename, size_t length, unsigned int patterns) {
    // If file already exists, skip
    if ( access(filename.c_str(), F_OK ) != -1 ) {
        std::cout << "Random patterns: " << filename << " already exist. Skip generating...\n";
        return;
    }
    std::ofstream fout(filename);
    int blockshift = 4;
    int twopower = 1 << blockshift;
    int maxRand = (1 << twopower) - 1;
    // use better pseudo-random generator than rand() % 2
    std::mt19937 e(std::random_device{}());
    std::uniform_int_distribution<int> uniform_dist(0, maxRand);
    // generate enought bits to p
    size_t numBlock = ((length - 1) >> blockshift) + 1;
    std::string p(numBlock << blockshift, '0');

    std::cout << "Generating " << patterns << " random patterns in file: " << filename << " ...\n";

    for (unsigned int i = 0; i < patterns; ++i) {
        for (size_t j = 0; j < numBlock; ++j) {
            int r = uniform_dist(e);
            // extract 8 bits from r
            for (int k = 0; k < twopower; ++k) {
                p[(j << blockshift) + k] = (r & 1) ? '1' : '0';
                r >>= 1;
            }
        }

        fout << p.substr(0, length) << std::endl;
    }
    fout.close();
}

void dumpLowNodes(Graph *G, std::string nodefilename)
{
    std::ofstream fnodes(nodefilename.c_str());
    for(auto &edg : G->lowprobEdges) {
        fnodes << edg->name << "\t" << edg->lowProbVal << "\t" << min(edg->signal0Prob, edg->signal1Prob) << std::endl;
    }
    fnodes.close();
}

void readLowNodes(Graph *G, std::string nodefilename)
{
    std::ifstream fnodes(nodefilename.c_str());
    while(fnodes.good()) {
        std::string name;
        fnodes >> name;
        if (name.length() == 0) break;
        Edge *edg = G->findEdge(name);
        if (edg == nullptr)
            cerr << "Edge " << name << "not found\n";
        fnodes >> edg->lowProbVal;
        if (edg->lowProbVal == 1)
            fnodes >> edg->signal1Prob;
        else if (edg->lowProbVal == 0)
            fnodes >> edg->signal0Prob;
        else
        {
            cerr << "error in signal val\n";
        }
        
        G->lowprobEdges.push_back(edg);
    }
    fnodes.close();
}
