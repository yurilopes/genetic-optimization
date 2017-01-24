#include <list>
#include <vector>
#include <bitset>
#include <random>
#include <ctime>

#include "population.h"

mt19937 gen0_1(static_cast<unsigned int>(std::time(0)));
uniform_real_distribution<float> disR0_1(0.0, 1.0);

#define randomReal() disR0_1(gen0_1)
#define randomBinary() 0

class GeneticAlgorithm {
	//Population *population;

	public:
		static void selectionRoulette(Population &pop);
		static Population* generateMatingPool(Population &pop);

	protected:
		static Individual * getFirstIndividualMatingRoulette(Population &pop, float probability);
		static Individual * getSecondIndividualMatingRoulette(Population &pop, Individual * ind0, float probability);
		

};


void GeneticAlgorithm::selectionRoulette(Population &pop) {
	int32_t total = 0, lowest, current;
	float accumulatedNormFit = 0.0;

	vector<Individual*>* population = pop.getIndividualVector();

	/*
	Calculate the lowest fitness value
	This is used to translate every fitness value to a positive range
	*/
	for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Individual *individual = *itr;
		current = individual->getFitness();
		if (itr == population->begin()) //Lowest = fitness of the first iteration
			lowest = current;
		if (current < lowest)
			lowest = current;
	}

	/*
	1 is added to the lowest value so that the worst individual can actually have a nonzero chance of mating
	*/
	lowest *= -1;
	lowest++;

	/*
	Redefine the fitness for each individual
	and calculate the total sum of fitness
	*/
	for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Individual *individual = *itr;
		total += individual->addToFitness(lowest);
	}

	/*
	Calculate the normalized fitness for each individual
	*/
	for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Individual *individual = *itr;
		accumulatedNormFit += ((float)individual->getFitness() / (float)total);		
		individual->setAccNormalizedFitness(accumulatedNormFit);
	}

}

Population* GeneticAlgorithm::generateMatingPool(Population & pop)
{
	Population *matingPopulation = new Population(0, 0); //Seeding values aren't important	

	int populationSize = pop.getIndividualVector()->size();

	for (int i = 0; i < populationSize; i += 2) {
		Individual *ind0 = getFirstIndividualMatingRoulette(pop, randomReal());
		Individual *ind1 = getSecondIndividualMatingRoulette(pop, ind0, randomReal());
		matingPopulation->getIndividualVector()->push_back(ind0);
		matingPopulation->getIndividualVector()->push_back(ind1);
	}

	return matingPopulation;
}

inline Individual * GeneticAlgorithm::getFirstIndividualMatingRoulette(Population & pop, float probability)
{
	vector<Individual*>* population = pop.getIndividualVector();

	for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Individual *individual = *itr;
		if (individual->getAccNormalizedFitness() >= probability)
			return individual;			
	}

	return pop.getFittestIndividual();
}

inline Individual * GeneticAlgorithm::getSecondIndividualMatingRoulette(Population & pop, Individual * ind0, float probability)
{
	vector<Individual*>* population = pop.getIndividualVector();

	for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Individual *individual = *itr;
		if ((individual->getAccNormalizedFitness() >= probability) && (individual != pop.getFittestIndividual()))
			return individual;
	}

	if ((ind0 != pop.getFittestIndividual()))
		return pop.getFittestIndividual();
	else
		return (*(std::next(population->begin()))); //Second element in list
}



void crossOver(vector<int32_t> &parentX, vector<int32_t> &parentY){
	/*
	Uniform crossover
	More info at: https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)
	*/
    int32_t vectorSize = parentX.size();
    vector<int32_t> *childX = new vector<int32_t>(vectorSize);
    vector<int32_t> *childY = new vector<int32_t>(vectorSize);

    bool first = false;

    for(int32_t i = 0; i<vectorSize; i++){
        bitset<32> bitX(parentX[i]);
        bitset<32> bitY(parentY[i]);

        for(int j=0; j<32; j++){
            /*
            The first bit is always kept untouched
            This is a premise of crossover
            */
            if ((i==0) && (j==0))
                continue;

            if(randomBinary()){
                /*
                Set bitX[j] as bitY[j] and vice-versa
                */
                bool bit = bitX.test(j);
                bitX.set(j, bitY.test(j));
                bitY.set(j, bit);
            }

        }

        int32_t x = (int32_t)bitX.to_ulong();
        int32_t y = (int32_t)bitY.to_ulong();

        (*childX)[i]=x;
        (*childY)[i]=y;

    }

}
