//Comment this line below if not running in Visual Studio
#include "stdafx.h"

using namespace std; 

#include "genetic-algorithm.h"

double fitnessFunction2(Chromosome * chromosome) {
	//http://tutorial.math.lamar.edu/Classes/Alg/NonlinearSystems.aspx

	float x = (*chromosome->getGenes())[0]->getValue().floatValue;
	float y = (*chromosome->getGenes())[1]->getValue().floatValue;

	//We have to check if the float/double values are valid
	if (isnan(x) || isinf(x))
		return -9999;
	if (isnan(y) || isinf(y))
		return -9999;


	double fitness = 10;
	
	if (2*x*x + y*y != 24)
		fitness -= abs(24 - (2 * x*x + y*y));
	
	if (x*x - y*y != -12)
		fitness -= abs(-12 - (x*x - y*y));

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

#define TOLERANCE			1e-4f
#define IDEAL_FITNESS		115.0f
#define MIN_SEED			-6.0f
#define MAX_SEED			6.0f
#define POPULATION_SIZE		8000
#define ITERATION_SHOW		1

int main(){

	/*
	First we create our chromosome model's genes

	We are declaring below that our chromosomes will look like this:

	CHROMOSOME = [UINT8 variable] [UINT8 variable] [UINT8 variable] [UINT8 variable]
	*/
	vector<Gene *> genes;

	Gene *gene = new Gene(FLOAT);
	gene->setSeedRange(float(MIN_SEED), float(MAX_SEED));
	genes.push_back(gene);
	
	gene = new Gene(FLOAT);
	gene->setSeedRange(float(MIN_SEED), float(MAX_SEED));
	genes.push_back(gene);

/*	gene = new Gene(FLOAT);
	gene->setSeedRange(float(0), float(MAX_SEED));
	genes.push_back(gene);

	gene = new Gene(FLOAT);
	gene->setSeedRange(float(0), float(MAX_SEED));
	genes.push_back(gene);
	*/
	/*
	Then we begin our GA setup
	*/
		
	GeneticAlgorithm ga(&genes);	

	ga.setElitism(true);	
	ga.setEliteSize(25);
	
	ga.setMutation(true);
	ga.setMutationRate(0.1f);
	
	ga.setFitnessFunction(fitnessFunction2);

	ga.initializePopulation(POPULATION_SIZE);

	clock_t timeBegin = clock(); //Starting time

	uint64_t i;
	for (i = 0; i < 0xFFFFFFFF; i++) {						

		ga.calculateFitness();		

		if (i % ITERATION_SHOW == 0) {
			cout << "Iteration " << i << endl;
			cout << "Fittest individual:" << endl;
			ga.printFittestChromosome();
		}
		//system("pause");
		if (abs(ga.getFittestChromosome()->getFitness() - IDEAL_FITNESS) <= TOLERANCE)
			break;

		ga.selectionRoulette();	
		ga.generateRouletteMatingPool();
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
