#include <list>
#include <vector>
#include <bitset>
#include <random>
#include <ctime>
#include <limits>

#include "population.h"

mt19937 genGA(static_cast<unsigned int>(std::time(0)));
uniform_real_distribution<float> disR0_1(0.0, 1.0);
uniform_real_distribution<float> disR0_100(0.0, 100.0);
uniform_int_distribution<int32_t> disI0_1(0, 1);

#define randomReal() disR0_1(genGA)
#define randomReal100() disR0_100(genGA)
#define randomBinary() disI0_1(genGA)

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
		void	printFittestIndividual();
		void	refreshFitnessFunction();
		void	setElitism(bool elit);
		bool	getElitism();
		bool	getMutation();
		void	setMutation(bool mut);
		float	getMutationRate();
		void	setMutationRate(float rate);
		void	setEliteSize(uint32_t amount);
		uint32_t	getEliteSize();
		void	printMatingPool();
		Chromosome * getFittestIndividual();


		void selectionRoulette();
		void generateRouletteMatingPool();
		void crossOver();

		static int32_t getVariable(vector<int32_t> *variables, uint32_t position);

		~GeneticAlgorithm();		

	protected:
		Population * gPopulation = NULL, * gMatingPool = NULL;
		int32_t minSeed = LONG_MIN, maxSeed = LONG_MAX;
		uint32_t populationSize = 0 , individualSize = 0;
		FitnessFunction fitnessFunc = NULL;		
		bool elitism = true;
		uint32_t eliteSize = 1;
		bool mutation = true;
		float mutationRate = 0.01f;
		
		void crossOverIndividualInt32(Individual * parentX, Individual * parentY, Individual * childXP, Individual * childYP);
		static Chromosome * getFirstIndividualMatingRoulette(Population * pop, float probability);
		static Chromosome * getSecondIndividualMatingRoulette(Population * pop, Individual * ind0, float probability);		
		void mutate(Individual * individual);		
};

inline void GeneticAlgorithm::selectionRoulette() {
	int32_t total = 0, lowest, current;
	float accumulatedNormFit = 0.0;

	vector<Chromosome*>* population = gPopulation->getIndividualVector();

	/*
	Calculate the lowest fitness value
	This is used to translate every fitness value to a positive range
	*/
	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *individual = *itr;
		current = individual->getFitness();
		if (itr == population->begin()) //Lowest = fitness of the first iteration
			lowest = current;
		if (current < lowest)
			lowest = current;
	}

	/*
	1 is added to the lowest value so that the worst individual can actually have a nonzero chance of mating
	*/
	if(lowest < 0)
		lowest *= -1;
	lowest++;

	/*
	Redefine the fitness for each individual
	and calculate the total sum of fitness
	*/
	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *individual = *itr;
		total += individual->addToAccFitness(lowest);
	}

	/*
	Calculate the normalized fitness for each individual
	*/
	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *individual = *itr;
		accumulatedNormFit += ((float)individual->getAccFitness() / (float)total);		
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

	if (eliteSize > popSize)
		return;

	vector<Chromosome*>* population = gPopulation->getIndividualVector();	
	
	if (elitism) {
		uint32_t i = 0;
		
		for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {			
			Chromosome *individual = *itr;
			Chromosome * clone = new Chromosome(individual);
			if (i >= eliteSize)
				break;
			gMatingPool->getIndividualVector()->push_back(clone);							
			i++;
		}
		size -= eliteSize;
	}
	
	for (uint32_t i = 0; i < size; i += 2) {
		Chromosome * ind0;
		Chromosome * clone0;						

		if (gMatingPool->getIndividualVector()->size() >= popSize) {
			break;
		}
		
		if (i+1 >= size) {
			//This individual will mate with the previous one
			//So we have to make sure it is NOT equal to the previous one
			ind0 = gMatingPool->getIndividualVector()->back(); //ind0 comes from Mating Pool
			ind0 = gMatingPool->getEqualIndividual(ind0); //ind0 comes from Population
			Chromosome * ind1 = getSecondIndividualMatingRoulette(gPopulation, ind0, randomReal());
			Chromosome * clone1 = new Chromosome(ind1);
			gMatingPool->getIndividualVector()->push_back(clone1);
			break;

		} else {			
			ind0 = getFirstIndividualMatingRoulette(gPopulation, randomReal());
			clone0 = new Chromosome(ind0);
			gMatingPool->getIndividualVector()->push_back(clone0);
		}
	
		Chromosome * ind1 = getSecondIndividualMatingRoulette(gPopulation, ind0, randomReal());
		Chromosome * clone1 = new Chromosome(ind1);
		gMatingPool->getIndividualVector()->push_back(clone1);
	}	

}

inline void GeneticAlgorithm::crossOver()
{
	if (gMatingPool == NULL)
		return;

	uint32_t popSize = gMatingPool->getIndividualVector()->size();	
	uint32_t startIndex = 0;
	
	if (elitism) {
		if (eliteSize > popSize)
			return; //This should never happen		
		startIndex = eliteSize;		
	}

	for (uint32_t i = startIndex; i < popSize; i += 2) {		
		//Iterate over the elements two at a time for crossover operation

		Chromosome * parentX = (*gMatingPool->getIndividualVector())[i];
		Chromosome * parentY;

		Chromosome * childXP = (*gPopulation->getIndividualVector())[i];
		Chromosome * childYP;

		if (i + 1 >= popSize) {
			parentY = (*gMatingPool->getIndividualVector())[i - 1]; //This guy is alone, so it'll mate with the previous Individual
			childYP = (*gPopulation->getIndividualVector())[i - 1];
		}
		else {
			parentY = (*gMatingPool->getIndividualVector())[i + 1];
			childYP = (*gPopulation->getIndividualVector())[i + 1];			
		}				

		//TODO Treat the situation when i+1 >= popsize takes over the i-1 already constructed child
		//This isn't really and issue, maybe this isn't a problem at all

		//Crossover each int32 in each Individual
		crossOverIndividualInt32(parentX, parentY, childXP, childYP);		
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

inline void GeneticAlgorithm::printFittestIndividual()
{
	if (gPopulation == NULL)
		return;
	gPopulation->getFittestIndividual()->print();
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

inline bool GeneticAlgorithm::getMutation()
{
	return mutation;
}

inline void GeneticAlgorithm::setMutation(bool mut)
{
	mutation = mut;
}

inline float GeneticAlgorithm::getMutationRate()
{
	return mutationRate;
}

inline void GeneticAlgorithm::setMutationRate(float rate)
{
	mutationRate = rate;
}

inline void GeneticAlgorithm::setEliteSize(uint32_t amount)
{
	eliteSize = amount;
}

inline uint32_t GeneticAlgorithm::getEliteSize()
{
	return eliteSize;
}

inline void GeneticAlgorithm::printMatingPool()
{
	if (gMatingPool != NULL)
		gMatingPool->printPopulation();
}

inline Chromosome * GeneticAlgorithm::getFittestIndividual()
{
	if (gPopulation == NULL)
		return NULL;
	return gPopulation->getFittestIndividual();
}

inline Chromosome * GeneticAlgorithm::getFirstIndividualMatingRoulette(Population * pop, float probability)
{
	vector<Chromosome*>* population = pop->getIndividualVector();

	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *individual = *itr;
		if (individual->getAccNormalizedFitness() >= probability)
			return individual;			
	}

	return population->back(); //This should only happen if the last individual has a normalized fitness less than 1
}

inline Chromosome * GeneticAlgorithm::getSecondIndividualMatingRoulette(Population * pop, Chromosome * ind0, float probability)
{

	vector<Chromosome*>* population = pop->getIndividualVector();

	bool found = false;

	Chromosome * result = NULL, * lastIndividual = NULL;

	//Generate a population that consists of the original excluding ind0
	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *individual = *itr;
		if (!individual->equals(ind0)) {
			lastIndividual = individual;			
			found = true;

			if (individual->getAccNormalizedFitness() >= probability && result == NULL) //Result was found
				result = individual;			
		}
	}

	if (!found) {
		//No individual different from ind0 was found
		//Population is already homogeneous		
		return ind0;
	}	

	if (result)
		return result;

	//tmpPop has individuals which are not ind0
	//The probability could not match any individuals in tmpPop
	//This happens if, and only if, the original population had its tail full of ind0s
	//In this case, we can either get the last individual in tmpPop which has the closest AccNormalizedFitness to probability
	//Or choose a random individual in tmpPop
	//The former, however, is empirically preferred since it proportionally favors the best individuals

	return lastIndividual;
}

inline void GeneticAlgorithm::mutate(Chromosome * individual)
{	
	if (mutationRate < randomReal100())
		return;

	uint32_t vectorSize = individual->getGenes()->size();

	uniform_int_distribution<uint32_t> dis2(0, 20);
	
	for (uint32_t i = 0; i < vectorSize; i++) {
		(*individual->getGenes())[i] = dis2(genGA);
	}
}

inline int32_t GeneticAlgorithm::getVariable(vector<int32_t>* variables, uint32_t position)
{
	return (*variables)[position];
}



void GeneticAlgorithm::crossOverIndividualInt32(Chromosome * parentX, Chromosome * parentY, Chromosome * childXP, Chromosome * childYP){
	/*
	Uniform crossover
	More info at: https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)

	This function changes both parents to the result of the crossover operation
	*/
	vector<int32_t> * parentXV = parentX->getGenes();
	vector<int32_t> * parentYV = parentY->getGenes();
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
	childXP->setChromosome(childX);
	if(mutation)
		mutate(childXP);
	childYP->setChromosome(childY);
	if(mutation)
		mutate(childYP);

}
