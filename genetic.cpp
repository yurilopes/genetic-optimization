//Comment this line below if not running in Visual Studio
#include "stdafx.h"

using namespace std; 

#include "genetic-algorithm.h"
#include "mutation-uniform.h"
#include "mutation-gaussian.h"
#include "mutation-vector.h"
#include "crossover-uniform-bitwise.h"

#include "examples.h"

#define ERROR_FITNESS		1e-4

#define	OPT_MODE			MODE_MAXIMIZE
#define OPTIMAL_FITNESS		-4.579582
#define POPULATION_SIZE		2000
#define ITERATION_SHOW		100
#define ELITE_SIZE			25
#define CROSSOVER_PROB		1.0f
#define MUTATION_PROB		0.05f
#define ES_NOFFSPRING		20
#define ES_ELITEONLY		true

int main(){	

	GeneticAlgorithm ga(getGenotype5());
	ga.setFitnessFunction(fitnessFunction5);
	ga.setOptimizationMode(OPT_MODE);
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
	MutationGaussian mutG(0.0, 0.01);
	MutationUniform mutU;
	vector<MutationVectorized *> mutvec;	
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	MutationVector mut(mutvec);
	

	ga.setMutationOperator(&mutU);
	ga.enableMutation(true);
	ga.setMutationProbability(MUTATION_PROB);	

	ga.initializePopulation(POPULATION_SIZE);

	clock_t timeBegin = clock(); //Starting time

	uint64_t i;
	bool eliteOnly = ES_ELITEONLY;
	for (i = 0; i < 0xFFFFFFFF; i++) {			

		ga.calculateFitness();			


		if (i % ITERATION_SHOW == 0) {
			cout << "Iteration " << i << endl;
			cout << "Fittest chromosome:" << endl;
			ga.printFittestChromosome();
		}

		if (abs(ga.getFittestChromosome()->getFitness() - OPTIMAL_FITNESS) <= ERROR_FITNESS)
			break;

		if (i < 2000) {
			ga.selectionRoulette();
			ga.generateRouletteMatingPool();
			ga.crossOver();
		}

		ga.evolutionStrategy(ES_NOFFSPRING, eliteOnly);

		if (i == 500)
			mutG.setStdDev(0.0001);		
		if (i == 1300) 
			mutG.setStdDev(0.00001);
		if (i == 2000) {
			mutG.setStdDev(0.000001);
			eliteOnly = false;
		}
		if (i == 4000) 
			mutG.setStdDev(0.0000001);					
		
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
