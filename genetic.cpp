//Comment this line below if not running in Visual Studio
#include "stdafx.h"

#include <limits.h>

using namespace std; 

//#include "genetic-algorithm.h"
#include "chromosome.h"


#define MAX_SEED 20

double fitnessFunction(Chromosome * chromosome){
	//http://www.zweigmedia.com/RealWorld/simplex.html

	/*int32_t x = GeneticAlgorithm::getVariable(variables, 0);
	int32_t y = GeneticAlgorithm::getVariable(variables, 1);
	int32_t z = GeneticAlgorithm::getVariable(variables, 2);
	int32_t w = GeneticAlgorithm::getVariable(variables, 3);

	int32_t fitness = (x / 2) + 3 * y + z + 4 * w;
	//Optimal Solution: p = 115; x = 10, y = 10, z = 0, w = 20

	bool violated = false;
	int32_t violations[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	/*
	To evaluate the constraints we check if their complement is true
	If so, the constraint has been violated and we should punish the fitness
	*/

/*	if (x + y + z + w > 40) {
		violated = true;
		violations[0] = (x + y + z + w) - 40;
	}
	if (2 * x + y - z - w < 10) {
		violated = true;
		violations[1] = 10 - (2 * x + y - z - w);
	}
	if (w - y < 10) {
		violated = true;
		violations[2] = 10 - (w - y);
	}

	if (x < 0) {
		violated = true;
		violations[3] = -x;
	}
	if (y < 0) {
		violated = true;
		violations[4] = -y;
	}
	if (z < 0) {
		violated = true;
		violations[5] = -z;
	}
	if (w < 0) {
		violated = true;
		violations[6] = -w;
	}
	if (x > MAX_SEED) {
		violated = true;
		violations[7] = x-MAX_SEED;
	}
	if (y > MAX_SEED) {
		violated = true;
		violations[8] = y - MAX_SEED;
	}
	if (z > MAX_SEED) {
		violated = true;
		violations[9] = z - MAX_SEED;
	}
	if (w > MAX_SEED) {
		violated = true;
		violations[10] = w - MAX_SEED;
	}

	if (violated) {		
		fitness = 0;
		for (int i = 0; i < 7; i++)
			fitness -= abs(violations[i]);
		for (int i = 3; i < 11; i++)
			if (violations[i] != 0)
				return -9999;
	}

	return fitness;
	*/
	return 0;
}

#define IDEAL_FITNESS 14

int main(){
    /*
    A population is vector of individuals
    Each individual is vector of ints
    Each problem variable/gene is an int
    */

	vector<Gene *> genes;

	Gene *gene = new Gene(UINT8);
	gene->setSeedRange(uint8_t(0), uint8_t(20));	
	genes.push_back(gene);
	
	gene = new Gene(INT8);
	gene->setSeedRange(int8_t(-50), int8_t(50));
	genes.push_back(gene);
	
	
	gene = new Gene(FLOAT);
	gene->setSeedRange(float(0), float(21));
	
	genes.push_back(gene);
	
	
	Chromosome ch(&genes, fitnessFunction);	

	for (size_t i = 0; i < genes.size(); i++) {
		Gene *g = (*ch.getGenes())[i];		
		if (g->getDataType() == FLOAT)
			cout << i << ": " << g->getDataType() << " {" << g->getValueFloat() << "} - [" << dec << (int16_t)g->getMinimumSeed().int8Value << ", " << g->getMaximumSeed().floatValue;					
		else
			cout << i << ": " << g->getDataType() << " {" << (int16_t)g->getValueInt8() << "} - [" << dec << (int16_t)g->getMinimumSeed().int8Value << ", " << (int16_t)g->getMaximumSeed().int8Value;
		cout << "]" << endl;
		
	}

	return 0;
/*
	GeneticAlgorithm ga;

	ga.setElitism(true);	
	ga.setEliteSize(5);
	ga.setMutation(true);
	ga.setMutationRate(0.01f);
	ga.setMinSeed(0);
	ga.setMaxSeed(MAX_SEED);
	//ga.setFitnessFunction(fitnessFunction);

	ga.initializePopulation(50, 4);

	clock_t timeBegin = clock(); //Starting time

	int i;
	for (i = 0; i < 0x6FFFFFFF; i++) {				

		if (i % 50 == 0) {
			cout << "Iteration " << i << endl;
			cout << "Fittest individual:" << endl;
			ga.printFittestIndividual();					
		}

		ga.calculateFitness();

		if (ga.getFittestIndividual()->getFitness() == IDEAL_FITNESS) 
			break;

		ga.selectionRoulette();	

		//cout << "Iteration " << i << endl;
		//ga.printPopulation();		

		//cout << "-------------------------------------" << endl;

		ga.generateRouletteMatingPool();
		//ga.printMatingPool();

		//cout << endl << endl;

		//system("pause");
		ga.crossOver();
	}

	clock_t timeEnd = clock(); //Ending time
	double timeSpent = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

	cout << "-------------------------------" << endl << endl;
	cout << endl << "Finished at iteration " << i << endl;
	cout << "Elapsed time: " << timeSpent << endl;
	cout << endl << "Final fittest individual: " << endl;
	ga.calculateFitness();
	ga.printFittestIndividual();
	cout << endl;
	//ga.printPopulation();

	

    return ga.getFittestIndividual()->getFitness();
	*/
}
