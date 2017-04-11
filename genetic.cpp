//Comment this line below if not running in Visual Studio
#include "stdafx.h"

using namespace std; 

#include "genetic-algorithm.h"
#include "mutation-uniform.h"
#include "mutation-gaussian.h"
#include "mutation-vector.h"
#include "crossover-uniform-bitwise.h"

#include "examples.h"

#define ERROR_FITNESS		1e-4f

#define OPTIMAL_FITNESS		-1.07654f
#define POPULATION_SIZE		2000
#define ITERATION_SHOW		10
#define ELITE_SIZE			25
#define CROSSOVER_PROB		1.0f
#define MUTATION_PROB		0.05f
#define ES_NOFFSPRING		4
#define ES_ELITEONLY		false

int main(){	

	GeneticAlgorithm ga(getGenotype3());
	ga.setFitnessFunction(fitnessFunction3);

	ga.setElitism(true);	
	ga.setEliteSize(ELITE_SIZE);
	
	CrossoverUniformBitwise cross;
	ga.setCrossoverOperator(&cross);
	ga.setCrossoverProbability(CROSSOVER_PROB);
		
	/*
	MutationGaussian mutGauss(0.0, 0.005);
	MutationUniform mutUniform;
	vector<MutationVectorized *> mutVect;
	mutVect.push_back(&mutGauss);
	mutVect.push_back(&mutUniform);
	MutationVector mut(mutVect); 
	*/

	//MutationUniform mut;
	MutationGaussian mut(0.0, 0.0001);
	
	ga.setMutationOperator(&mut);
	ga.enableMutation(true);
	ga.setMutationProbability(MUTATION_PROB);	

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

		if (abs(ga.getFittestChromosome()->getFitness() - OPTIMAL_FITNESS) <= ERROR_FITNESS)
			break;

		ga.selectionRoulette();	
		ga.generateRouletteMatingPool();
		ga.crossOver();

		ga.evolutionStrategy(ES_NOFFSPRING, ES_ELITEONLY);

		if (i == 400)
			mut.setStdDev(0.00001);
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

	

    return (int)ga.getFittestChromosome()->getFitness();
	
}
