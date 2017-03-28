//Comment this line below if not running in Visual Studio
#include "stdafx.h"

using namespace std; 

#include "genetic-algorithm.h"
#include "mutation-uniform.h"
#include "mutation-gaussian.h"
#include "crossover-uniform-bitwise.h"

#define ERROR_FITNESS			1e-5f
#define BOOST_FACTOR			1e3f
#define	PUNISHMENT_FACTOR		1e3f

double fitnessFunction2(Chromosome * chromosome) {
	//http://www.zweigmedia.com/RealWorld/simplex.html

	float x = (*chromosome->getGenes())[0]->getValue().floatValue;
	uint8_t y = (*chromosome->getGenes())[1]->getValue().uint8Value;
	//float z = (*chromosome->getGenes())[2]->getValue().floatValue;
	//float w = (*chromosome->getGenes())[3]->getValue().floatValue;

	//We have to check if the float/double values are valid
	if (isnan(x) || isinf(x))
		return -INFINITY;
/*	if (isnan(y) || isinf(y))
		return -9999;
/*	if (isnan(z) || isinf(z))
		return -9999;
	if (isnan(w) || isinf(w))
		return -9999;
*/

	double fitness = -BOOST_FACTOR*(2 * x + y);	

	bool violated = false;
	float violations[4] = { 0, 0, 0, 0 };

	/*
	For constraint evaluation, we check if their complement is true.
	If so, the constraint has been violated and we should punish the fitness proportionately.
	*/

	if (1.25f > x*x + (float)y) {
		violated = true;
		violations[0] = 1.25f - (x*x + (float)y);
	}
	if (x + y > 1.6f) {
		violated = true;
		violations[1] = 1.6f - (x+y);
	}
	if (x > 1.6f) {
		violated = true;
		violations[2] = 1.6f - x;
	}

	if (x < 0) {
		violated = true;
		violations[3] = -x;
	}

	if (violated) {
		fitness = 0;
		for (int i = 0; i <4; i++)
			fitness -= abs(violations[i]);
		fitness *= PUNISHMENT_FACTOR;
	}

	return fitness;
}

double fitnessFunction(Chromosome * chromosome){	
	//http://www.zweigmedia.com/RealWorld/simplex.html

	float x = (*chromosome->getGenes())[0]->getValue().floatValue;
	float y = (*chromosome->getGenes())[1]->getValue().floatValue;
	float z = (*chromosome->getGenes())[2]->getValue().floatValue;
	float w = (*chromosome->getGenes())[3]->getValue().floatValue;
	
	//We have to check if the float/double values are valid
	if (isnan(x) || isinf(x))
		return -9999;
	if (isnan(y) || isinf(y))
		return -9999;
	if (isnan(z) || isinf(z))
		return -9999;
	if (isnan(w) || isinf(w))
		return -9999;
	

	double fitness = (x / 2) + 3 * y + z + 4 * w;
	//Optimal Solution: p = 115; x = 10, y = 10, z = 0, w = 20

	bool violated = false;
	float violations[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	/*
	For constraint evaluation, we check if their complement is true.
	If so, the constraint has been violated and we should punish the fitness proportionately.
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

#define OPTIMAL_FITNESS		2.0f
#define MIN_SEED			0.0f
#define MAX_SEED			1.6f
#define POPULATION_SIZE		2000
#define ITERATION_SHOW		1

int main(){	

	/*
	First we create our genotype

	We are declaring below that our chromosomes will look like this:

	CHROMOSOME = [UINT8 variable] [UINT8 variable] [UINT8 variable] [UINT8 variable]
	*/
	vector<Gene *> genotype;

	/*
	Gene *gene = new Gene(UINT8);
	gene->setSeedRange(0, 40);
	genotype.push_back(gene);
	*/

	GeneValue minSeed, maxSeed, minSeed2, maxSeed2;
	minSeed.floatValue = MIN_SEED;
	maxSeed.floatValue = MAX_SEED;

	maxSeed2.uint8Value = 1;
	minSeed2.uint8Value = 1;
	
	Gene *gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	genotype.push_back(gene);

	gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype.push_back(gene);

	/*gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	genotype.push_back(gene);

	gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	genotype.push_back(gene);
	*/
	
	/*
	Then we begin our GA setup
	*/
		
	GeneticAlgorithm ga(&genotype);

	ga.setElitism(true);	
	ga.setEliteSize(25);

	
	CrossoverUniformBitwise cross;
	ga.setCrossoverOperator(&cross);
	ga.setCrossoverProbability(1.0f);
	
	//vector<char> mutEnabled = { true, true, false, true };
	//MutationGaussian mut(0.0, 0.005, &mutEnabled);
	MutationGaussian mut(0.0, 0.005);
	ga.setMutationOperator(&mut);
	ga.enableMutation(true);
	ga.setMutationProbability(0.01f);	
	
	ga.setFitnessFunction(fitnessFunction2);

	ga.initializePopulation(POPULATION_SIZE);

	clock_t timeBegin = clock(); //Starting time

	uint64_t i;
	for (i = 0; i < 0xFFFFFFFF; i++) {			

		ga.calculateFitness();	

		if (i % ITERATION_SHOW == 0) {
			cout << "Iteration " << i << endl;
			cout << "Fittest chromosome:" << endl;
			ga.printFittestChromosome();
		}
		//system("pause");
		if (abs(ga.getFittestChromosome()->getFitness() - OPTIMAL_FITNESS) <= ERROR_FITNESS)
			break;

		ga.selectionRoulette();	
		ga.generateRouletteMatingPool();
		ga.crossOver();
	}

	clock_t timeEnd = clock(); //Ending time
	double timeSpent = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

	cout << "-------------------------------" << endl << endl;
	cout << endl << "Finished at iteration " << i << endl;
	cout << "Elapsed time: " << timeSpent << "s" << endl;
	cout << endl << "Final fittest chromosome: " << endl;
	ga.calculateFitness();
	ga.printFittestChromosome();
	cout << endl;
	//ga.printPopulation();

	

    return (int)ga.getFittestChromosome()->getFitness();
	
}
