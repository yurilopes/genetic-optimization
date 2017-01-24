#include <list>
#include <vector>
#include <bitset>
#include <random>
#include <ctime>
#include <limits>
#include <Windows.h>

#include "population.h"

mt19937 gen0_1(static_cast<unsigned int>(std::time(0)));
uniform_real_distribution<float> disR0_1(0.0, 1.0);
uniform_int_distribution<int32_t> disI0_1(0, 1);

#define randomReal() disR0_1(gen0_1)
#define randomBinary() disI0_1(gen0_1)

class GeneticAlgorithm {
	public:
		int32_t	getMinSeed();
		void	setMinSeed(int32_t mSeed);
		int32_t	getMaxSeed();
		void	setMaxSeed(int32_t mSeed);
		uint32_t getIndividualSize();
		uint32_t getPopulationSize();		
		void	setFitnessFunction(FitnessFunction fitFunc);
		FitnessFunction getFitnessFunction();
		void	initializePopulation(uint32_t populationSize, uint32_t individualSize);
		void	calculateFitness();
		void	printPopulation();		
		void	refreshFitnessFunction();
		void	setElitism(bool elit);
		bool	getElitism();
		void	setEliteAmount(uint32_t amount);
		uint32_t	getEliteAmount();
		void	printMatingPool();

		void selectionRoulette();
		void generateRouletteMatingPool();
		void crossOver();

		~GeneticAlgorithm();

	protected:
		Population * gPopulation = NULL, * gMatingPool = NULL;
		int32_t minSeed = LONG_MIN, maxSeed = LONG_MAX;
		uint32_t populationSize = 0 , individualSize = 0;
		FitnessFunction fitnessFunc = NULL;		
		bool elitism = true;
		uint32_t eliteAmout = 1;
		
		static void crossOverIndividualInt32(Individual * parentX, Individual * parentY);
		static Individual * getFirstIndividualMatingRoulette(Population * pop, float probability);
		static Individual * getSecondIndividualMatingRoulette(Population * pop, Individual * ind0, float probability);
};

inline void GeneticAlgorithm::selectionRoulette() {
	int32_t total = 0, lowest, current;
	float accumulatedNormFit = 0.0;

	vector<Individual*>* population = gPopulation->getIndividualVector();

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

void GeneticAlgorithm::generateRouletteMatingPool()
{
	Sleep(100);
	if (gMatingPool != NULL)
		delete gMatingPool;
		
	gMatingPool = new Population(0, 5); //Seeding values aren't important
	gMatingPool->setMainPopulation(false); //The destructor won't try to free invalid Individual pointers

	unsigned int popSize = gPopulation->getIndividualVector()->size();
	uint32_t size = popSize;

	if (eliteAmout > popSize)
		return;


	vector<Individual*>* population = gPopulation->getIndividualVector();	
	
	if (elitism) {
		uint32_t i = 0;
		
		for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {			
			Individual *individual = *itr;
			if (i >= eliteAmout)
				break;
			gMatingPool->getIndividualVector()->push_back(individual);							
			i++;
		}
		size -= eliteAmout;
	}
	
	for (uint32_t i = 0; i < size; i += 2) {
		Individual *ind0 = getFirstIndividualMatingRoulette(gPopulation, randomReal());
		Individual *ind1 = getSecondIndividualMatingRoulette(gPopulation, ind0, randomReal());		
		gMatingPool->getIndividualVector()->push_back(ind0);		
		if (gMatingPool->getIndividualVector()->size() >= popSize)
			break;
		gMatingPool->getIndividualVector()->push_back(ind1);
	}	

}

inline void GeneticAlgorithm::crossOver()
{
	if (gMatingPool == NULL)
		return;

	uint32_t popSize = gMatingPool->getIndividualVector()->size();	
	uint32_t startIndex = 0;
	
	if (elitism) {
		if (eliteAmout > popSize)
			return; //This should never happen
		startIndex = eliteAmout;		
	}
	



	for (uint32_t i = startIndex; i < popSize; i += 2) {		
		//Iterate over the elements two at a time for crossover operation
		int k;
		Individual * parentX = (*gMatingPool->getIndividualVector())[i];
		Individual * parentY;
		if (i + 1 >= popSize) {
			parentY = (*gMatingPool->getIndividualVector())[i - 1]; //This guy is alone, so it'll mate with the previous Individual
			k = i - 1;
		}
		else {
			parentY = (*gMatingPool->getIndividualVector())[i + 1];
			k = i + 1;
		}

		cout << i << " & " << k << endl;
		//cin >> k;

		//Crossover each int32 in each Individual
		crossOverIndividualInt32(parentX, parentY);
		//printMatingPool();
	}
	
}

inline GeneticAlgorithm::~GeneticAlgorithm()
{
	if (gPopulation != NULL)
		delete gPopulation;
	if (gMatingPool != NULL)
		delete gMatingPool;	
}


inline int32_t GeneticAlgorithm::getMinSeed()
{
	return minSeed;
}

inline void GeneticAlgorithm::setMinSeed(int32_t mSeed)
{
	minSeed = mSeed;
}

inline int32_t GeneticAlgorithm::getMaxSeed()
{
	return maxSeed;
}

inline uint32_t GeneticAlgorithm::getIndividualSize()
{
	return individualSize;
}

inline uint32_t GeneticAlgorithm::getPopulationSize()
{
	return populationSize;
}

inline void GeneticAlgorithm::setFitnessFunction(FitnessFunction fitFunc)
{
	fitnessFunc = fitFunc;
}

inline FitnessFunction GeneticAlgorithm::getFitnessFunction()
{
	return fitnessFunc;
}

inline void GeneticAlgorithm::initializePopulation(uint32_t populationSize, uint32_t individualSize)
{
	if (gPopulation != NULL)
		delete gPopulation;

	gPopulation = new Population(minSeed, maxSeed);

	if (fitnessFunc != NULL)
		gPopulation->initialize(populationSize, individualSize, fitnessFunc);
}

inline void GeneticAlgorithm::setMaxSeed(int32_t mSeed)
{
	maxSeed = mSeed;
}

inline void GeneticAlgorithm::calculateFitness()
{
	if (gPopulation != NULL)
		gPopulation->calculateFitness();
}

inline void GeneticAlgorithm::printPopulation()
{
	if (gPopulation != NULL)
		gPopulation->printPopulation();
}

inline void GeneticAlgorithm::refreshFitnessFunction()
{
	if (gPopulation == NULL)
		return;

	gPopulation->refreshFitnessFunction(fitnessFunc);
}

inline void GeneticAlgorithm::setElitism(bool elit)
{
	elitism = elit;
}

inline bool GeneticAlgorithm::getElitism()
{
	return elitism;
}

inline void GeneticAlgorithm::setEliteAmount(uint32_t amount)
{
	eliteAmout = amount;
}

inline uint32_t GeneticAlgorithm::getEliteAmount()
{
	return eliteAmout;
}

inline void GeneticAlgorithm::printMatingPool()
{
	if (gMatingPool != NULL)
		gMatingPool->printPopulation();
}

inline Individual * GeneticAlgorithm::getFirstIndividualMatingRoulette(Population * pop, float probability)
{
	vector<Individual*>* population = pop->getIndividualVector();

	for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Individual *individual = *itr;
		if (individual->getAccNormalizedFitness() >= probability)
			return individual;			
	}

	return pop->getFittestIndividual();
}

inline Individual * GeneticAlgorithm::getSecondIndividualMatingRoulette(Population * pop, Individual * ind0, float probability)
{
	vector<Individual*>* population = pop->getIndividualVector();

	for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Individual *individual = *itr;
		if ((individual->getAccNormalizedFitness() >= probability) && (individual != pop->getFittestIndividual()))
			return individual;
	}

	if ((ind0 != pop->getFittestIndividual()))
		return pop->getFittestIndividual();
	else
		return (*(std::next(population->begin()))); //Second element in list
}



void GeneticAlgorithm::crossOverIndividualInt32(Individual * parentX, Individual * parentY){
	/*
	Uniform crossover
	More info at: https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)

	This function changes both parents to the result of the crossover operation
	*/
	vector<int32_t> * parentXV = parentX->getIndividual();
	vector<int32_t> * parentYV = parentY->getIndividual();
    int32_t vectorSize = parentXV->size();
    vector<int32_t> *childX = new vector<int32_t>(vectorSize);
    vector<int32_t> *childY = new vector<int32_t>(vectorSize);

    bool first = false;

    for(int32_t i = 0; i<vectorSize; i++){
        bitset<32> bitX((*parentXV)[i]);
        bitset<32> bitY((*parentYV)[i]);

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

	//Assign the childs to replace their parents
	parentX->setIndividual(childX);
	parentY->setIndividual(childY);

}
