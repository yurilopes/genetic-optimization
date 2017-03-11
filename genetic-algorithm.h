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
		void				setFitnessFunction(FitnessFunction fitFunc);
		FitnessFunction		getFitnessFunction();
		void				initializePopulation(uint32_t populationSize);
		void				calculateFitness();
		void				printPopulation();		
		void				printFittestChromosome();
		void				refreshFitnessFunction();
		void				setElitism(bool elit);
		bool				getElitism();
		bool				getMutation();
		void				setMutation(bool mut);
		float				getMutationRate();
		void				setMutationRate(float rate);
		void				setEliteSize(uint32_t amount);
		uint32_t			getEliteSize();
		size_t				getPopulationSize();
		void				printMatingPool();
		Chromosome			*getFittestChromosome();


		void				selectionRoulette();
		void				generateRouletteMatingPool();
		void				crossOver();

		static int32_t getVariable(vector<int32_t> *variables, uint32_t position);

		GeneticAlgorithm(vector<Gene *>	*genModel);
		~GeneticAlgorithm();

	protected:
		vector<Gene *>		*geneModel;
		Population			*gPopulation = NULL, * gMatingPool = NULL;		
		uint32_t			populationSize = 0;
		FitnessFunction		fitnessFunc = NULL;		
		bool				elitism = true;
		uint32_t			eliteSize = 1;
		bool				mutation = true;
		float				mutationRate = 0.01f;
		
		void crossOverUniform(Chromosome * parentX, Chromosome * parentY, Chromosome * childXP, Chromosome * childYP);
		static Chromosome * getFirstChromosomeMatingRoulette(Population * pop, float probability);
		static Chromosome * getSecondChromosomeMatingRoulette(Population * pop, Chromosome * ind0, float probability);		
		void mutateUniform(Chromosome * Chromosome);		
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
	if(lowest < 0)
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
		
	gMatingPool = new Population(geneModel);
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
		
		if (i+1 >= size) {			
			//This Chromosome will mate with the previous one
			//So we have to make sure it is NOT equal to the previous one
			ind0 = gMatingPool->getChromosomes()->back(); //ind0 comes from Mating Pool			
			ind0 = gMatingPool->getEqualChromosome(ind0); //ind0 comes from Population
			Chromosome * ind1 = getSecondChromosomeMatingRoulette(gPopulation, ind0, randomReal());
			Chromosome * clone1 = new Chromosome(ind1);
			gMatingPool->getChromosomes()->push_back(clone1);
			break;
		} else {			
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

		//Crossover each int32 in each Chromosome
		crossOverUniform(parentX, parentY, childXP, childYP);		
	}
	
}

inline GeneticAlgorithm::GeneticAlgorithm(vector<Gene*>* genModel)
{
	geneModel = genModel;
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

	gPopulation = new Population(geneModel);	
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

inline Chromosome * GeneticAlgorithm::getFittestChromosome()
{
	if (gPopulation == NULL)
		return NULL;
	return gPopulation->getFittestChromosome();
}

inline Chromosome * GeneticAlgorithm::getFirstChromosomeMatingRoulette(Population * pop, float probability)
{
	vector<Chromosome*>* population = pop->getChromosomes();

	for (vector<Chromosome*>::iterator itr = population->begin(); itr != population->end(); itr++) {
		Chromosome *Chromosome = *itr;
		if (Chromosome->getAccNormalizedFitness() >= probability)
			return Chromosome;			
	}

	return population->back(); //This should only happen if the last Chromosome has a normalized fitness less than 1
}

inline Chromosome * GeneticAlgorithm::getSecondChromosomeMatingRoulette(Population * pop, Chromosome * ind0, float probability)
{

	vector<Chromosome*>* population = pop->getChromosomes();

	bool found = false;

	Chromosome * result = NULL, * lastChromosome = NULL;

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

inline void GeneticAlgorithm::mutateUniform(Chromosome * chromosome)
{	
	if (mutationRate < randomReal100())
		return;

	vector<Gene *> *genes = chromosome->getGenes();
	for (vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {
		Gene * gene = *it;

		switch (gene->getDataType())
		{
			case (FLOAT): {
				uniform_real_distribution<float> randomF(gene->getMinimumSeed().floatValue, gene->getMaximumSeed().floatValue);
				gene->setValueFloat(randomF(genGA));
				break;
			}
			case (DOUBLE): {
				uniform_real_distribution<double> randomD(gene->getMinimumSeed().doubleValue, gene->getMaximumSeed().doubleValue);
				gene->setValueDouble(randomD(genGA));
				break;
			}
			case (INT8): {
				uniform_int_distribution<int16_t> randomI16(gene->getMinimumSeed().int8Value, gene->getMaximumSeed().int8Value);
				#pragma warning(push)
				#pragma warning(disable : 4244)
				/*
				Should never cause loss of data, since seed values are within the 8 bit range
				*/
				gene->setValueUInt8(randomI16(genGA));
				#pragma warning(pop)
				break;
			}
			case (UINT8): {
				/*
				Can't have a uniform_int_distribution<uint8_t> because:

				It seems that this library implementaion behaviour does not violate the ISO C++ standard, because (26.5.1.1)
				"effect of instantiating a template ... that has a template type parameter named IntType __is undefined unless__ the
				corresponding template argument is cv-unqualified
				and __is one of short, int, long, long long, unsigned short, unsigned int, unsigned long, or unsigned long long__",
				char is not listed.

				More info at: http://stackoverflow.com/questions/31460733/why-arent-stduniform-int-distributionuint8-t-and-stduniform-int-distri
				*/
				uniform_int_distribution<uint16_t> randomI16(gene->getMinimumSeed().uint8Value, gene->getMaximumSeed().uint8Value);
				#pragma warning(push)
				#pragma warning(disable : 4244)
				/*
				Should never cause loss of data, since seed values are within the 8 bit range
				*/

				gene->setValueUInt8(randomI16(genGA));
				#pragma warning(pop)
				break;
			}
			case (INT16): {
				uniform_int_distribution<int16_t> randomI16(gene->getMinimumSeed().int16Value, gene->getMaximumSeed().int16Value);
				gene->setValueInt16(randomI16(genGA));
				break;
			}
			case (UINT16): {
				uniform_int_distribution<uint16_t> randomI16(gene->getMinimumSeed().uint16Value, gene->getMaximumSeed().uint16Value);
				gene->setValueUInt16(randomI16(genGA));
				break;
			}
			case (INT32): {
				uniform_int_distribution<int32_t> randomI32(gene->getMinimumSeed().int32Value, gene->getMaximumSeed().int32Value);
				gene->setValueInt32(randomI32(genGA));
				break;
			}
			case (UINT32): {
				uniform_int_distribution<uint32_t> randomI32(gene->getMinimumSeed().uint32Value, gene->getMaximumSeed().uint32Value);
				gene->setValueUInt32(randomI32(genGA));
				break;
			}
			case (INT64): {
				uniform_int_distribution<int64_t> randomI64(gene->getMinimumSeed().int64Value, gene->getMaximumSeed().int64Value);
				gene->setValueInt64(randomI64(genGA));
				break;
			}
			case (UINT64): {
				uniform_int_distribution<uint64_t> randomI64(gene->getMinimumSeed().uint64Value, gene->getMaximumSeed().uint64Value);
				gene->setValueUInt64(randomI64(genGA));
				break;
			}
			default: {
				//TODO: Set value for CUSTOM data type
				break;
			}
		}
	}
	
}

inline int32_t GeneticAlgorithm::getVariable(vector<int32_t>* variables, uint32_t position)
{
	return (*variables)[position];
}



void GeneticAlgorithm::crossOverUniform(Chromosome * parentX, Chromosome * parentY, Chromosome * childXP, Chromosome * childYP){
	/*
	Uniform crossover
	More info at: https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)

	This function changes both parents to the result of the crossover operation
	*/

	//First we clone the parents
	Chromosome *childX = new Chromosome(parentX);
	Chromosome *childY = new Chromosome(parentY);
	vector<Gene *> *childXGenes = childX->getGenes();
	vector<Gene *> *childYGenes = childY->getGenes();

    bool first = false;
	
	for (unsigned int i = 0; i < childXGenes->size(); i++) {
		//First we iterate through the genes

		Gene *genXV = (*childXGenes)[i];
		Gene *genYV = (*childYGenes)[i];

		for (unsigned int j = 0; j < genXV->getGeneBits()->size(); j++) {
			//Then we iterate through the bits of each gene pair
			if (!first) {
				/*
				The first bit is always kept untouched
				This is a premise of crossover
				*/
				first = true;
				continue;
			}

			if (randomBinary()) {
				/*
				Set bitX[j] as bitY[j] and vice-versa
				*/				
				unsigned char bitX = (*genXV->getGeneBits())[j];
				(*genXV->getGeneBits())[j] = (*genYV->getGeneBits())[j];
				(*genYV->getGeneBits())[j] = bitX;

			}

		}

	
	}

	//Assign the childs to replace their parents
	childXP->setChromosome(childX);
	if(mutation)
		mutateUniform(childXP);
	childYP->setChromosome(childY);
	if(mutation)
		mutateUniform(childYP);		
}
