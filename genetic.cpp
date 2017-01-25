//Comment this line below if running in Visual Studio
#include "stdafx.h"

#include <limits.h>

using namespace std; 

#include "genetic-algorithm.h"

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

int32_t fitnessFunction(vector<int32_t> *variables) {
	int32_t x0 = GeneticAlgorithm::getVariable(variables, 0);
	int32_t x1 = GeneticAlgorithm::getVariable(variables, 1);

	int32_t fitness = -2*x0 + 5*x1;	

	/*
	To evaluate the constraints we look if their complement is true
	If so, the constraint has been violated and we should punish the fitness
	*/
	if (x0 < 100)
		fitness = -9999; //Constraint violation, punish the fitness
	if (x0 > 200)
		fitness = -9999;
	if (x1 < 80)
		fitness = -9999;
	if (x1 > 170)
		fitness = -9999;
	if (x1 < -x0 + 200)
		fitness = -9999;

	return fitness;		
}

#define IDEAL_FITNESS 650

int main(){
    /*
    A population is vector of individuals
    Each individual is vector of ints
    Each problem variable/gene is an int
    */
	GeneticAlgorithm ga;

	ga.setElitism(true);	
	ga.setMutation(true);
	ga.setMutationRate(0.01f);
	ga.setMinSeed(0);
	ga.setMaxSeed(500);
	ga.setFitnessFunction(fitnessFunction);

	ga.initializePopulation(400, 2);

	clock_t timeBegin = clock(); //Starting time

	int i;
	for (i = 0; i < 1000; i++) {		
		if (i % 50 == 0 && i !=0) {
			cout << "Iteration " << i << endl;
			cout << "Fittest individual:" << endl;
			ga.printFittestIndividual();			
		}

		ga.calculateFitness();

		if (ga.getFittestIndividual()->getFitness() == IDEAL_FITNESS) 
			break;

		if (i == 0) {
			cout << "First fittest individual:" << endl;
			ga.printFittestIndividual();
			cout << "-------------------------------" << endl << endl;
		}

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
	ga.printFittestIndividual();	


	

    return ga.getFittestIndividual()->getFitness();
}
