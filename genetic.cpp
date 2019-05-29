/*
 * Authors: Yangdi Lyu
 * Jan 2018
 */

#include <algorithm>
#include "graph.h"
#include "genetic.h"
#include "traverse.h"

void geneOneVectorGA(Graph* G, string& origvec, string& pairvec) {
    // Start genetic algorithm
    srand(1234);

    Population pop(200, 5, origvec.length(), 0.1, origvec);
    string bestvec;

    for (int i = 0; i < pop.generation; ++i) {
        pop.getFitness(G, origvec);
        pop.selectPop();
        pop.crossover();
        pop.mutation();
    }

    pairvec = pop.bestvec;
}

Population::Population(int size, int gen, int len, double mrate, string& vec):
        popsize(size), generation(gen), indvlength(len), 
         bestfit(-1), mutrate(mrate) 
{
    indv.resize(size, vec);
    fitness.resize(size, -1);
    int diff = (vec.length() + 125) / 250;
    if (diff == 0)
        diff = 1;
    // Randomly flip bits for individuals
    for (int i = 0; i < size; ++i) {
        for (int j = 0 ; j < diff; ++j) {
            int idx = rand() % len;
            indv[i][idx] = (vec[idx] == '0') ? '1' : '0';
        }
    }
}

void Population::getFitness(Graph* G, string& vec)
{
    int i, singlefit, bestidx = -1;
    for (i = 0; i < popsize; i++) { 
        singlefit = (int) (10000 * switchSimPatternOnePair(G, indv[i], vec));
        if (singlefit > bestfit) {
            bestfit = singlefit;
            bestidx = i;
        }
        fitness[i] = singlefit;
    }
    // Cumulative fitnss
    for (i = 1; i < popsize; i++) { 
        fitness[i] += fitness[i - 1];
    }
    // If best finess change, change bestvec
    if (bestidx != -1) {
        bestvec = indv[bestidx];
    }
}

void Population::selectPop() {
    vector<string> newindv;
    int maxfit = fitness.back();
    for (int i = 0; i < popsize; ++i) {
        int r = rand() % (maxfit + 1);
        // Select based on fitness
        vector<int>::iterator it = lower_bound(fitness.begin(), fitness.end(), r);
        newindv.push_back(indv[it - fitness.begin()]);
    }
    indv = newindv;
}

void Population::crossover() {
    // Crossover each pair
    for (int i = 0; i + 1 < popsize; i += 2) {
        // Cross point
        int r = rand() % indvlength;
        for (int j = 0; j < r; j++) {
            char temp = indv[i][j];
            indv[i][j] = indv[i + 1][j]; 
            indv[i + 1][j] = temp; 
        }
    }
}

void Population::mutation() {
    for (int i = 0; i < popsize; i++) {
        int r = rand() % 100;
        if (r < mutrate * 100) {
            // Flip the j-th bit
            int j = rand() % indvlength;
            indv[i][j] = (indv[i][j] == '0')? '1' : '0';
        }
    }
}


