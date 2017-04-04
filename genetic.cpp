//Comment this line below if not running in Visual Studio
#include "stdafx.h"

using namespace std; 

#include "genetic-algorithm.h"
#include "mutation-uniform.h"
#include "mutation-gaussian.h"
#include "crossover-uniform-bitwise.h"

#include "examples.h"

#define ERROR_FITNESS		1e-3f

#define OPTIMAL_FITNESS		-2.124f
#define POPULATION_SIZE		4000
#define ITERATION_SHOW		10
#define ELITE_SIZE			50
#define CROSSOVER_PROB		1.0f
#define MUTATION_PROB		0.05f

int main(){	

	GeneticAlgorithm ga(getGenotype2d());
	ga.setFitnessFunction(fitnessFunction2d);

	ga.setElitism(true);	
	ga.setEliteSize(ELITE_SIZE);
	
	CrossoverUniformBitwise cross;
	ga.setCrossoverOperator(&cross);
	ga.setCrossoverProbability(CROSSOVER_PROB);
	
	//vector<char> mutEnabled = { true, true, false, true };
	//MutationGaussian mut(0.0, 0.005, &mutEnabled);
	MutationGaussian mut(0.0, 0.5);
	//MutationUniform mut;
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
