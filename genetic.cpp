//Comment this line below if not running in Visual Studio
#include "stdafx.h"

#include <limits.h>

using namespace std; 

#include "genetic-algorithm.h"

/*Na mating pool, garantir que um indivíduo nao vai cruzar com ele mesmo não só pelo ponteiro,
mas usando equals. se o laço principal nao der certo, tentar de novo com uma populacao clonada sem aquele individuo
se isso nao escolher ninguem, que ele cruze consigo porque a populacao ja esta condenada
*/
/*
Maximize F = –2*x[0] + 5*x[1]
Subject to:

x[0] >= 100
x[0] <= 200

x[1] >= 80
x[1] <= 170

x[1] >= –x[0] + 200

Optimal values:
P = 650
x[0] = 100
x[1] = 170

First example of http://www.purplemath.com/modules/linprog3.htm

*/
/*
int32_t fitnessFunction(vector<int32_t> *variables) {
	int32_t x0 = GeneticAlgorithm::getVariable(variables, 0);
	int32_t x1 = GeneticAlgorithm::getVariable(variables, 1);

	int32_t fitness = -2*x0 + 5*x1;	
	
	int violations = 0;

	/*
	To evaluate the constraints we look if their complement is true
	If so, the constraint has been violated and we should punish the fitness
	*/ /*
	if (x0 < 100)
		violations++;
	if (x0 > 200)
		violations++;
	if (x1 < 80)
		violations++;
	if (x1 > 170)
		violations++;
	if (x1 < -x0 + 200)
		violations++;

	if (violations > 0) {
		fitness = 0;
		for (int i = 0; i < violations; i++)
			fitness -= 5;
	}

	return fitness;		
}
*/

#define MAX_SEED 20


int32_t fitnessFunction(vector<int32_t> *variables) {
	//http://www.zweigmedia.com/RealWorld/simplex.html

	int32_t x = GeneticAlgorithm::getVariable(variables, 0);
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
}


#define IDEAL_FITNESS 115

int main(){
    /*
    A population is vector of individuals
    Each individual is vector of ints
    Each problem variable/gene is an int
    */
	GeneticAlgorithm ga;

	ga.setElitism(true);	
	ga.setEliteSize(5);
	ga.setMutation(true);
	ga.setMutationRate(0.01f);
	ga.setMinSeed(0);
	ga.setMaxSeed(MAX_SEED);
	ga.setFitnessFunction(fitnessFunction);

	ga.initializePopulation(2000, 4);

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
}
