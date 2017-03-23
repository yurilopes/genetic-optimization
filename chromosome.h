#pragma once

#include <vector>
#include <inttypes.h>
#include <ctime>
#include <random>
#include "gene.h"

#define MAXFLOATPRINT 999.f
#define MINFLOATPRINT -99.f

extern mt19937 genGA;

class Chromosome; //Forward declaration of Chromosome class for FitnessFunction typedef

typedef double(*FitnessFunction)(Chromosome * chromosome);

class Chromosome{
    protected:
        vector<Gene *>						*genes = NULL;
        double								fitness;
		double								accFitness;
		double								accNormalizedFitness;
        FitnessFunction						fitnessFunction;
		vector<Gene *>						*genotype;

		void								deleteGenes();
		

    public:
        Chromosome(vector<Gene *> * geneModel, FitnessFunction fitFunction);
		Chromosome(Chromosome * original);
        ~Chromosome();
        vector<Gene *>						*getGenes();
		vector<Gene *>						*getGenotype();
		double								calculateFitness();
		double								addToAccFitness(double value);
        double								getFitness();
		double								getAccFitness();
		void								setFitnessFunction(FitnessFunction fitFunc);
		FitnessFunction						getFitnessFunction();
		double								getAccNormalizedFitness();
		void								setAccNormalizedFitness(double fit);
		void								print();
		bool								equals(Chromosome * ind);

        static bool							compare(Chromosome *ind0, Chromosome *ind1);
};


inline void Chromosome::deleteGenes()
{
	if (!genes)
		return;
	for (vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {
		Gene * gene = *it;
		delete gene;
	}
	delete genes;
}

Chromosome::Chromosome(vector<Gene *> * genModel, FitnessFunction fitFunction){
	genotype = genModel;	
	deleteGenes();

	genes = new vector<Gene *>();
	
	for (unsigned int i = 0; i < genModel->size(); i++) {
		Gene * gModel = (*genModel)[i];
		Gene *newGene = new Gene(gModel);		
		newGene->initialize();
		//Add new gene to gene vector
		genes->push_back(newGene);
	}
    
    fitness = 0;
	accNormalizedFitness = 0.0;
    fitnessFunction = fitFunction;
}

inline Chromosome::Chromosome(Chromosome * original)
{
	fitness = original->getFitness();
	accNormalizedFitness = original->getAccNormalizedFitness();
	fitnessFunction = original->getFitnessFunction();

	/*
	Now we clone the gene vector.
	A simple
		genes = new vector<Gene *>(*original->getGenes());
	won't do!

	We have to create copies of each gene. Ohterwise we are just copying pointers and the objects themselves still have only one instance
		in memory.
	We want to create copies of the Gene objects' instances.
	*/

	deleteGenes();

	genes = new vector<Gene *>();
	for (vector<Gene *>::iterator it = original->getGenes()->begin(); it != original->getGenes()->end(); it++) {
		Gene *gen = *it;
		Gene *newGene = new Gene(gen);
		genes->push_back(newGene);
	}	
}

Chromosome::~Chromosome(){
	deleteGenes();
}

vector<Gene *> * Chromosome::getGenes(){
    return genes;	
}

inline vector<Gene*> * Chromosome::getGenotype()
{
	return genotype;
}

double Chromosome::getFitness(){
    return fitness;
}

inline double Chromosome::getAccFitness()
{
	return accFitness;
}

inline void Chromosome::setFitnessFunction(FitnessFunction fitFunc)
{
	fitnessFunction = fitFunc;
}

inline FitnessFunction Chromosome::getFitnessFunction()
{
	return fitnessFunction;
}

double Chromosome::getAccNormalizedFitness() {
	return accNormalizedFitness;
}

void Chromosome::setAccNormalizedFitness(double fit)
{
	accNormalizedFitness = fit;
}

inline void Chromosome::print()
{
	for (vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {
		Gene * gen = *it;
		switch (gen->getDataType()) {
			case INT8:	
				printf("%6" PRIi8 , gen->getValue().int8Value);
				break;
			case UINT8:
				printf("%6" PRIu8, gen->getValue().uint8Value);
				break;
			case INT16:
				printf("%6" PRIi16, gen->getValue().int16Value);
				break;
			case UINT16:
				printf("%6" PRIu16, gen->getValue().uint16Value);
				break;
			case INT32:
				printf("%6" PRIi32, gen->getValue().int32Value);
			case UINT32:
				printf("%6" PRIu32, gen->getValue().uint32Value);
				break;
			case INT64:
				printf("%6" PRIi64, gen->getValue().int64Value);
			case UINT64:
				printf("%6" PRIu64, gen->getValue().uint64Value);
				break;
			case FLOAT:
				if (gen->getValue().floatValue >= MINFLOATPRINT && gen->getValue().floatValue <= MAXFLOATPRINT)
					printf("%7.5f", gen->getValue().floatValue);
				else
					printf("%7.5e", gen->getValue().floatValue);
				break;
			case DOUBLE:		
				if (gen->getValue().doubleValue >= MINFLOATPRINT && gen->getValue().doubleValue <= MAXFLOATPRINT)
					printf("%7.5f", gen->getValue().doubleValue);
				else
					printf("%7.5e", gen->getValue().doubleValue);				
				break;		
		}		
		std::cout << ", \t";
	}
	if(fitness >= MINFLOATPRINT && fitness <= MAXFLOATPRINT)
		printf("F=%7.5f\n", fitness);
	else
		printf("F=%7.5e\n", fitness);
}

inline bool Chromosome::equals(Chromosome * ind)
{
	if (genes->size() != ind->getGenes()->size())		
		return false;	
	
	for (unsigned int i = 0; i < genes->size(); i++) {
		/*
		Compare the highest possible value
		It doesn't matter the data type, be it float, int8, int64 ou double
		This comparison should work given the way the data is stored in memory
		*/	
		if ((*genes)[i]->getDataType() != (*ind->getGenes())[i]->getDataType())
			return false;

		//Either double or uint64_t will be the biggest data field
		if (sizeof(double) > sizeof(uint64_t)) {
			if ((*genes)[i]->getValue().doubleValue != (*ind->getGenes())[i]->getValue().doubleValue)
				return false;
		}
		else {
			if ((*genes)[i]->getValue().uint64Value != (*ind->getGenes())[i]->getValue().uint64Value)
				return false;
		}						
	}
	return true;

}

double Chromosome::calculateFitness(){
	fitness = fitnessFunction(this);	
	accFitness = fitness;
    return fitness;
}

inline double Chromosome::addToAccFitness(double value)
{
	return accFitness+=value;
}

bool Chromosome::compare(Chromosome *ind0, Chromosome *ind1){
    /*
    Maximization of the fitness function is implied here
    */
    return ind0->getFitness() > ind1->getFitness();
}

