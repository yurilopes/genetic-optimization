//Comment this line below if not running in Visual Studio
//#include "stdafx.h"

using namespace std;

#include "genetic-algorithm.h"
#include "mutation-uniform.h"
#include "mutation-gaussian.h"
#include "mutation-vector.h"
#include "crossover-uniform-bitwise.h"

#include "examples.h"

#define FITNESS_CLOSENESS	0.9999

#define	OPT_MODE			MODE_MAXIMIZE
#define OPTIMAL_FITNESS		-285519.9508
#define POPULATION_SIZE		1000
#define ITERATION_SHOW		100
#define ELITE_SIZE			25
#define CROSSOVER_PROB		1.0f
#define MUTATION_PROB		0.05f
#define ES_NOFFSPRING		50
#define ES_ELITEONLY		false
#define IT_STOP_CROSSOVER	500

int main(){

	GeneticAlgorithm ga(getGenotype8());
	ga.setFitnessFunction(fitnessFunction8);
	ga.setOptimizationMode(OPT_MODE);
	ga.setElitism(true);
	ga.setEliteSize(ELITE_SIZE);

	CrossoverUniformBitwise cross;
	ga.setCrossoverOperator(&cross);
	ga.setCrossoverProbability(CROSSOVER_PROB);

	MutationGaussian mutG(0.0, 0.01);
	MutationUniform mutU;
	vector<MutationVectorized *> mutvec;
	mutvec.push_back(&mutU); //N
	mutvec.push_back(&mutU);
	mutvec.push_back(&mutU);
	mutvec.push_back(&mutU);
	mutvec.push_back(&mutU);
	mutvec.push_back(&mutU);
	mutvec.push_back(&mutG); //V
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG); //B
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG); //TL
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	mutvec.push_back(&mutG);
	MutationVector mut(mutvec);


	ga.setMutationOperator(&mutU);
	ga.enableMutation(true);
	ga.setMutationProbability(MUTATION_PROB);

	ga.initializePopulation(POPULATION_SIZE);

	clock_t timeBegin = clock(); //Starting time

	Chromosome *lastChm = new Chromosome(ga.getFittestChromosome());

	printf("CHAR_BIT %d\n", CHAR_BIT);

	uint64_t i, countEqual = 0;
	bool eliteOnly = ES_ELITEONLY;
	for (i = 0; i < 0xFFFFFFFF; i++) {

		ga.calculateFitness();


		if (i % ITERATION_SHOW == 0) {
			cout << "Iteration " << i << endl;
			cout << "Fittest chromosome:" << endl;
			ga.printFittestChromosome();
			cout << "Std Dev: " << mutG.getStdDev() << endl;
		}

		if (ga.getFittestChromosome()->getFitness() / OPTIMAL_FITNESS >= FITNESS_CLOSENESS
			&&
			ga.getFittestChromosome()->getFitness() / OPTIMAL_FITNESS <= 1.0 + (1.0 - FITNESS_CLOSENESS))
			break;

		//Last chromosome is not equals the best one
		if (!ga.getFittestChromosome()->equals(lastChm)) {
			delete lastChm;
			lastChm = new Chromosome(ga.getFittestChromosome());
			countEqual = 0;
		}
		else {
			countEqual++;
		}

		if (countEqual >= 200){
			countEqual = 0;
			if(mutG.getStdDev()>=0.001)
				mutG.setStdDev(mutG.getStdDev()/10.0);
		}

		if (i < IT_STOP_CROSSOVER) { //Stop crossover and let ES guide the population
			ga.selectionRoulette();
			ga.generateRouletteMatingPool();
			ga.setMutationOperator(&mutU);
			ga.crossOver();
		}

		ga.setMutationOperator(&mut);
		ga.evolutionStrategy(ES_NOFFSPRING, eliteOnly);

		if (i == IT_STOP_CROSSOVER) {
			ga.setPopulationSurvivalSize(ELITE_SIZE);
		}

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
