//Comment this line below if not running in Visual Studio
#include "stdafx.h"

#include <limits.h>

using namespace std; 

//#include "genetic-algorithm.h"
#include "genetic-algorithm.h"

double fitnessFunction(Chromosome * chromosome){	
	//http://www.zweigmedia.com/RealWorld/simplex.html

	uint8_t x = (*chromosome->getGenes())[0]->getValueUInt8();
	uint8_t y = (*chromosome->getGenes())[1]->getValueUInt8();
	uint8_t z = (*chromosome->getGenes())[2]->getValueUInt8();
	uint8_t w = (*chromosome->getGenes())[3]->getValueUInt8();
	

	double fitness = (x / 2) + 3 * y + z + 4 * w;
	//Optimal Solution: p = 115; x = 10, y = 10, z = 0, w = 20

	bool violated = false;
	int32_t violations[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	/*
	To evaluate the constraints we check if their complement is true
	If so, the constraint has been violated and we should punish the fitness
	*/

	if (x + y + z + w > 40) {		
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

	if (violated) {		
		fitness = 0;
		for (int i = 0; i < 7; i++)
			fitness -= abs(violations[i]);
	}

	return fitness;	
}

#define IDEAL_FITNESS 115
#define MAX_SEED 20
#define POPULATION_SIZE 200

int main(){

	/*
	First we create our chromosome model's genes

	We are declaring below that our chromosomes will look like this:

	CHROMOSOME = [UINT8 variable] [UINT8 variable] [UINT8 variable] [UINT8 variable]
	*/
	vector<Gene *> genes;

	Gene *gene = new Gene(UINT8);
	gene->setSeedRange(uint8_t(0), uint8_t(MAX_SEED));
	genes.push_back(gene);
	
	gene = new Gene(UINT8);
	gene->setSeedRange(uint8_t(0), uint8_t(MAX_SEED));
	genes.push_back(gene);

	gene = new Gene(UINT8);
	gene->setSeedRange(uint8_t(0), uint8_t(MAX_SEED));
	genes.push_back(gene);

	gene = new Gene(UINT8);
	gene->setSeedRange(uint8_t(0), uint8_t(MAX_SEED));
	genes.push_back(gene);

	/*
	Then we begin our GA setup
	*/
		
	GeneticAlgorithm ga(&genes);	

	ga.setElitism(true);	
	ga.setEliteSize(5);
	
	ga.setMutation(true);
	ga.setMutationRate(0.1f);
	
	ga.setFitnessFunction(fitnessFunction);

	ga.initializePopulation(POPULATION_SIZE);	

	clock_t timeBegin = clock(); //Starting time

	int i;
	for (i = 0; i < 0x6FFFFFFF; i++) {				

		if (i % 50 == 0) {
			cout << "Iteration " << i << endl;
			cout << "Fittest individual:" << endl;
			ga.printFittestChromosome();					
		}

		ga.calculateFitness();

		if (ga.getFittestChromosome()->getFitness() == IDEAL_FITNESS)
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
	ga.printFittestChromosome();
	cout << endl;
	//ga.printPopulation();

	

    return (int)ga.getFittestChromosome()->getFitness();
	
}
