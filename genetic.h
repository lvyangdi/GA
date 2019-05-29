/*
 * Authors: Yangdi Lyu
 * Jan 2018
 */

#ifndef __GENETIC_H__
#define __GENETIC_H__

#include <iostream>
#include <cstdlib>
#include <string>
#include <ctype.h>
#include <unistd.h>
using namespace std;

class Graph;

class Population{
public:
    int popsize;
    int generation;
    int indvlength;
    string bestvec;
    double bestfit;
    vector<string> indv;
    vector<int> fitness;   // cumulative fitness
    double mutrate;

    Population(int size, int gen, int len, double mrate, string& vec);
    void selectPop();
    void crossover();
    void mutation();
    void getFitness(Graph* G, string& vec);
};

void geneOneVectorGA(Graph* G, string& origvec, string& pairvec);

#endif
