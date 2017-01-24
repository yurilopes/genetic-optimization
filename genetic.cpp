//Comment this line below if running in Visual Studio
#include "stdafx.h"

#include <limits.h>

using namespace std; 

#include "genetic-algorithm.h"


int32_t fitnessFunction(vector<int32_t> *variables){
    int32_t soma = 0;
    for(vector<int32_t>::iterator it = variables->begin(); it!=variables->end(); it++){
        soma += *it;
    }
    return soma;
}

int main(){
    /*
    A population is list of individuals
    Each individual is vector of ints
    Each problem variable/gene is an int
    */
	GeneticAlgorithm ga;
	ga.setElitism(false);
	ga.setMinSeed(-5);
	ga.setMaxSeed(5);
	ga.setFitnessFunction(fitnessFunction);
	ga.initializePopulation(6, 3);
	ga.calculateFitness();			
	ga.selectionRoulette();
	ga.printPopulation();
	cout << "------------------" << endl;
	ga.generateRouletteMatingPool();
	ga.printMatingPool();

    return 0;
}
