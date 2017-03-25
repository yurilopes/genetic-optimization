#pragma once

#include <list>
#include <vector>
#include <bitset>
#include <random>
#include <ctime>
#include <limits>

#if CHAR_BIT != 8
#pragma message("CHAR_BIT is not 8 bits long. This could lead to erratic behaviour and loss of information.")
#endif // CHAR_BIT != 8

#include "population.h"
#include "mutation.h"
#include "crossover.h"

mt19937 genGA(static_cast<unsigned int>(std::time(0)));
uniform_real_distribution<float> disR0_1(0.0, 1.0);
uniform_int_distribution<int16_t> disI0_1(0, 1);

class GeneticAlgorithm {
public:
	void				setFitnessFunction(FitnessFunction fitFunc);
	FitnessFunction		getFitnessFunction();
	void				initializePopulation(uint32_t populationSize);
	void				calculateFitness();
	void				printPopulation();
	void				printFittestChromosome();
	void				refreshFitnessFunction();
	void				setElitism(bool elit);
	bool				getElitism();
	bool				isMutationEnabled();
	void				enableMutation(bool mut);
	float				getMutationProbability();
	void				setMutationProbability(float prob);
	double				getCrossoverProbability();
	void				setCrossoverProbability(float prob);
	void				setEliteSize(uint32_t amount);
	uint32_t			getEliteSize();
	size_t				getPopulationSize();
	void				printMatingPool();
	Chromosome			*getFittestChromosome();
	vector<Gene *>		*getGenotype();
	void				setMutationOperator(Mutation *mut);
	Mutation			*getMutationOperator();
	void				setCrossoverOperator(Crossover *cross);
	Crossover			*getCrossoverOperator();


	void				selectionRoulette();
	void				generateRouletteMatingPool();
	void				crossOver();

	static float		randomReal();
	static uint16_t		randomBinary();

	GeneticAlgorithm(vector<Gene *>	*gentyp);
	~GeneticAlgorithm();

protected:
	vector<Gene *>		*genotype;
	Population			*gPopulation = NULL, *gMatingPool = NULL;
	uint32_t			populationSize = 0;
	FitnessFunction		fitnessFunc = NULL;
	Crossover			*crossoverOperator = NULL;
	Mutation			*mutationOperator = NULL;
	bool				elitism = true;
	uint32_t			eliteSize = 1;
	bool				mutationEnabled = false;
	float				mutationProbability = 0.01f;
	float				crossoverProbability = 1.0f;

	static Chromosome * getFirstChromosomeMatingRoulette(Population * pop, double probability);
	static Chromosome * getSecondChromosomeMatingRoulette(Population * pop, Chromosome * ind0, double probability);
	void				mutateChromosome(Chromosome * Chromosome);
};

inline void GeneticAlgorithm::selectionRoulette() {
	double total = 0, lowest, current;
	double accumulatedNormFit = 0.0;

	vector<Chromosome*>* chromosomes = gPopulation->getChromosomes();

	/*
	Calculate the lowest fitness value
	This is used to translate every fitness value to a positive range
	*/
	for (vector<Chromosome*>::iterator itr = chromosomes->begin(); itr != chromosomes->end(); itr++) {
		Chromosome *Chromosome = *itr;
		current = Chromosome->getFitness();
		if (itr == chromosomes->begin()) //Lowest = fitness of the first iteration
			lowest = current;
		if (current < lowest)
			lowest = current;
	}

	/*
	1 is added to the lowest value so that the worst Chromosome can actually have a nonzero chance of mating
	*/
	if (lowest < 0)
		lowest *= -1;
	lowest++;

	/*
	Redefine the fitness for each Chromosome
	and calculate the total sum of fitness
	*/
	for (vector<Chromosome*>::iterator itr = chromosomes->begin(); itr != chromosomes->end(); itr++) {
		Chromosome *Chromosome = *itr;
		total += Chromosome->addToAccFitness(lowest);
	}

	/*
	Calculate the normalized fitness for each Chromosome
	*/
	for (vector<Chromosome*>::iterator itr = chromosomes->begin(); itr != chromosomes->end(); itr++) {
		Chromosome *Chromosome = *itr;
		accumulatedNormFit += (Chromosome->getAccFitness() / total);
		Chromosome->setAccNormalizedFitness(accumulatedNormFit);
	}

}

inline void GeneticAlgorithm::generateRouletteMatingPool()
{
	if (gMatingPool != NULL)
		delete gMatingPool;

	gMatingPool = new Population(genotype);
	gMatingPool->refreshFitnessFunction(fitnessFunc);

	size_t popSize = gPopulation->getChromosomes()->size();
	size_t size = popSize;

	if (eliteSize > popSize)
		return;

	vector<Chromosome*>* population = gPopulation->getChromosomes();

	if (elitism) {
		size_t i = 0;

		for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
			if (i >= eliteSize)
				break;
			Chromosome *chromosome = *itr;
			Chromosome * clone = new Chromosome(chromosome);
			gMatingPool->getChromosomes()->push_back(clone);
			i++;
		}
		size -= eliteSize;
	}

	for (size_t i = 0; i < size; i += 2) {
		Chromosome * ind0;
		Chromosome * clone0;

		if (gMatingPool->getChromosomes()->size() >= popSize) {
			break;
		}

		if (i + 1 >= size) {
			//This Chromosome will mate with the previous one
			//So we have to make sure it is NOT equal to the previous one
			ind0 = gMatingPool->getChromosomes()->back(); //ind0 comes from Mating Pool			
			ind0 = gMatingPool->getEqualChromosome(ind0); //ind0 comes from Population
			Chromosome * ind1 = getSecondChromosomeMatingRoulette(gPopulation, ind0, randomReal());
			Chromosome * clone1 = new Chromosome(ind1);
			gMatingPool->getChromosomes()->push_back(clone1);
			break;
		}
		else {
			ind0 = getFirstChromosomeMatingRoulette(gPopulation, randomReal());
			clone0 = new Chromosome(ind0);
			gMatingPool->getChromosomes()->push_back(clone0);
		}

		Chromosome * ind1 = getSecondChromosomeMatingRoulette(gPopulation, ind0, randomReal());
		Chromosome * clone1 = new Chromosome(ind1);
		gMatingPool->getChromosomes()->push_back(clone1);
	}

}

inline void GeneticAlgorithm::crossOver()
{
	if (gMatingPool == NULL)
		return;

	size_t popSize = gMatingPool->getChromosomes()->size();
	size_t startIndex = 0;

	if (elitism) {
		if (eliteSize > popSize)
			return; //This should never happen		
		startIndex = eliteSize;
	}

	for (size_t i = startIndex; i < popSize; i += 2) {
		//Iterate over the elements two at a time for crossover operation

		Chromosome * parentX = (*gMatingPool->getChromosomes())[i];
		Chromosome * parentY;

		Chromosome * childXP = (*gPopulation->getChromosomes())[i];
		Chromosome * childYP;

		if (i + 1 >= popSize) {
			parentY = (*gMatingPool->getChromosomes())[i - 1]; //This guy is alone, so it'll mate with the previous Chromosome
			childYP = (*gPopulation->getChromosomes())[i - 1];
		}
		else {
			parentY = (*gMatingPool->getChromosomes())[i + 1];
			childYP = (*gPopulation->getChromosomes())[i + 1];
		}

		//Crossover both parent chromosomes
		if (randomReal() <= crossoverProbability) {
			crossoverOperator->crossover(parentX, parentY, childXP, childYP);

			//If mutation is enabled, try to mutate the children
			if (mutationEnabled) {
				mutateChromosome(childXP);
				mutateChromosome(childYP);
			}
		}
	}

}

inline float GeneticAlgorithm::randomReal()
{
	return disR0_1(genGA);
}

inline uint16_t GeneticAlgorithm::randomBinary()
{
	return disI0_1(genGA);
}

inline GeneticAlgorithm::GeneticAlgorithm(vector<Gene*>* gentyp)
{
	genotype = gentyp;
}

inline GeneticAlgorithm::~GeneticAlgorithm()
{
	if (gPopulation != NULL)
		delete gPopulation;
	if (gMatingPool != NULL)
		delete gMatingPool;
}

inline size_t GeneticAlgorithm::getPopulationSize()
{
	if (gPopulation && gPopulation->getChromosomes())
		return gPopulation->getChromosomes()->size();
}

inline void GeneticAlgorithm::setFitnessFunction(FitnessFunction fitFunc)
{
	fitnessFunc = fitFunc;
	if (gPopulation)
		gPopulation->refreshFitnessFunction(fitFunc);
}

inline FitnessFunction GeneticAlgorithm::getFitnessFunction()
{
	return fitnessFunc;
}

inline void GeneticAlgorithm::initializePopulation(uint32_t populationSize)
{
	if (gPopulation != NULL)
		delete gPopulation;

	gPopulation = new Population(genotype);
	gPopulation->initialize(populationSize, fitnessFunc);
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

inline void GeneticAlgorithm::printFittestChromosome()
{
	if (gPopulation == NULL)
		return;
	gPopulation->getFittestChromosome()->print();
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

inline bool GeneticAlgorithm::isMutationEnabled()
{
	return mutationEnabled;
}

inline void GeneticAlgorithm::enableMutation(bool mut)
{
	mutationEnabled = mut;
}

inline float GeneticAlgorithm::getMutationProbability()
{
	return mutationProbability;
}

inline void GeneticAlgorithm::setMutationProbability(float prob)
{
	mutationProbability = prob;
}

inline double GeneticAlgorithm::getCrossoverProbability()
{
	return crossoverProbability;
}

inline void GeneticAlgorithm::setCrossoverProbability(float prob)
{
	crossoverProbability = prob;
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

inline Chromosome * GeneticAlgorithm::getFittestChromosome()
{
	if (gPopulation == NULL)
		return NULL;
	return gPopulation->getFittestChromosome();
}

inline vector<Gene*>* GeneticAlgorithm::getGenotype()
{
	return genotype;
}

inline void GeneticAlgorithm::setMutationOperator(Mutation * mut)
{
	mutationOperator = mut;
}

inline Mutation * GeneticAlgorithm::getMutationOperator()
{
	return mutationOperator;
}

inline void GeneticAlgorithm::setCrossoverOperator(Crossover * cross)
{
	crossoverOperator = cross;
}

inline Crossover * GeneticAlgorithm::getCrossoverOperator()
{
	return crossoverOperator;
}

inline Chromosome * GeneticAlgorithm::getFirstChromosomeMatingRoulette(Population * pop, double probability)
{
	vector<Chromosome*>* population = pop->getChromosomes();

	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *Chromosome = *itr;
		if (Chromosome->getAccNormalizedFitness() >= probability)
			return Chromosome;
	}

	return population->back(); //This should only happen if the last Chromosome has a normalized fitness less than 1
}

inline Chromosome * GeneticAlgorithm::getSecondChromosomeMatingRoulette(Population * pop, Chromosome * ind0, double probability)
{

	vector<Chromosome*>* population = pop->getChromosomes();

	bool found = false;

	Chromosome * result = NULL, *lastChromosome = NULL;

	//Generate a population that consists of the original excluding ind0
	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *Chromosome = *itr;
		if (!Chromosome->equals(ind0)) {
			lastChromosome = Chromosome;
			found = true;

			if (Chromosome->getAccNormalizedFitness() >= probability && result == NULL) //Result was found
				result = Chromosome;
		}
	}

	if (!found) {
		//No Chromosome different from ind0 was found
		//Population is already homogeneous		
		return ind0;
	}

	if (result)
		return result;

	//tmpPop has Chromosomes which are not ind0
	//The probability could not match any Chromosomes in tmpPop
	//This happens if, and only if, the original population had its tail full of ind0s
	//In this case, we can either get the last Chromosome in tmpPop which has the closest AccNormalizedFitness to probability
	//Or choose a random Chromosome in tmpPop
	//The former, however, is empirically preferred since it proportionally favors the best Chromosomes

	return lastChromosome;
}

inline void GeneticAlgorithm::mutateChromosome(Chromosome * chromosome)
{

	if (randomReal() <= mutationProbability)
		mutationOperator->mutate(chromosome);

}