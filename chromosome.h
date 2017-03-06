#pragma once

#include <vector>
#include <conio.h>
#include <inttypes.h>
#include "gene.h"

typedef double (*FitnessFunction)(Chromosome &chromosome);

class Chromosome{
    protected:
        vector<Gene>		*genes;
        double				fitness;
		double				accFitness;
		double				accNormalizedFitness;
        FitnessFunction		fitnessFunction;

    public:
        Chromosome(const vector<Gene> &geneModel, FitnessFunction fitFunction);
		Chromosome(Chromosome * original);
        ~Chromosome();
        vector<Gene>* getGenes();
		void setChromosome(Chromosome &chm);
		double calculateFitness();
		double addToAccFitness(double value);
        double getFitness();
		double getAccFitness();
		void setFitnessFunction(FitnessFunction fitFunc);
		FitnessFunction getFitnessFunction();
		double getAccNormalizedFitness();
		void setAccNormalizedFitness(double fit);
		void print();
		bool equals(Chromosome * ind);

        static bool compare(Chromosome *ind0, Chromosome *ind1);
};


Chromosome::Chromosome(const vector<Gene> &geneModel, FitnessFunction fitFunction){
	//TODO
    genes = new vector<Gene>(indSize);
    fitness = 0;
	accNormalizedFitness = 0.0;
    fitnessFunction = fitFunction;
}

inline Chromosome::Chromosome(Chromosome * original)
{
	fitness = original->getFitness();
	accNormalizedFitness = original->getAccNormalizedFitness();
	fitnessFunction = original->getFitnessFunction();

	//Clone the vector
	genes = new vector<Gene>(*original->getGenes());
}

Chromosome::~Chromosome(){
    delete genes;
}

vector<Gene>* Chromosome::getGenes(){
    return genes;
}

inline void Chromosome::setChromosome(Chromosome &chm)
{
	vector<Gene> * ind = chm.getGenes;
	if (genes != NULL)
		delete genes;
	genes = ind;
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
	for (vector<Gene>::iterator it = genes->begin(); it != genes->end(); it++) {
		Gene gen = *it;
		switch (gen.getDataType) {
			case INT8:				
				printf("%6" PRIi8 " (%6" PRIu8 ")", gen.getValueInt8, gen.getValueUInt8);
				break;
			case INT16:
				printf("%6" PRIi16 " (%6" PRIu16 ")", gen.getValueInt16, gen.getValueUInt16);
				break;
			case INT32:
				printf("%6" PRIi32 " (%6" PRIu32 ")", gen.getValueInt32, gen.getValueUInt32);
				break;
			case INT64:
				printf("%6" PRIi64 " (%6" PRIu64 ")", gen.getValueInt64, gen.getValueUInt64);
				break;
			case FLOAT:
				printf("%6.3f", gen.getValueFloat);
				break;
			case DOUBLE:			
				printf("%6.3f", gen.getValueDouble);
				break;
			case CUSTOM:
				//TODO
				break;			
		}		
		cout << ", \t";
	}
	cout << "F=" << fitness << endl;;	
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
		if ((*genes)[i].getValueUInt64 != (*ind->getGenes())[i].getValueUInt64)
			return false;
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

