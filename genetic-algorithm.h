#include <list>
#include <vector>
#include <bitset>
#include <random>
#include <ctime>
#include <limits>

#include "population.h"

mt19937 gen0_1(static_cast<unsigned int>(std::time(0)));
uniform_real_distribution<float> disR0_1(0.0, 1.0);

#define randomReal() disR0_1(gen0_1)
#define randomBinary() 0

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

		~GeneticAlgorithm();

	protected:
		Population * gPopulation = NULL, * gMatingPool = NULL;
		int32_t minSeed = LONG_MIN, maxSeed = LONG_MAX;
		uint32_t populationSize = 0 , individualSize = 0;
		FitnessFunction fitnessFunc = NULL;		
		bool elitism = true;
		uint32_t eliteAmout = 1;

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

inline void GeneticAlgorithm::generateRouletteMatingPool()
{
	if (gMatingPool != NULL)
		delete gMatingPool;

	gMatingPool = new Population(0, 0); //Seeding values aren't important

	unsigned int popSize = gPopulation->getIndividualVector()->size();
	uint32_t size = popSize;

	if (eliteAmout > popSize)
		return;


	//vector<Individual*>* population = gPopulation->getIndividualVector();	
	/*
	if (elitism) {
		uint32_t i = 0;
		for (vector<Individual*>::iterator itr = population->begin(); itr != population->end(); itr++) {
			Individual *individual = *itr;
			gMatingPool->getIndividualVector()->push_back(individual);
			cout << i << endl;
			i++;
			if (i >= eliteAmout)
				break;
		}
		size -= eliteAmout;
	}
	*/

	cout << "Elite ok" << endl;

	//return;



	for (uint32_t i = 0; i < size; i += 2) {
		Individual *ind0 = getFirstIndividualMatingRoulette(gPopulation, randomReal());
		Individual *ind1 = getSecondIndividualMatingRoulette(gPopulation, ind0, randomReal());
		gMatingPool->getIndividualVector()->push_back(ind0);
		gMatingPool->getIndividualVector()->push_back(ind1);
	}
	/*
	while (gMatingPool->getIndividualVector()->size() < popSize) {
		cout << "Here" << endl;
		Individual *ind0 = getFirstIndividualMatingRoulette(gPopulation, randomReal());
		gMatingPool->getIndividualVector()->push_back(ind0);
	}*/

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
