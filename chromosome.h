#pragma once

#include <vector>
#include <conio.h>
#include <inttypes.h>
#include <ctime>
#include <random>
#include "gene.h"
#include "chromosome.h"

class Chromosome; //Forward declaration of Chromosome class for FitnessFunction typedef

mt19937 genGA(static_cast<unsigned int>(std::time(0)));
typedef double(*FitnessFunction)(Chromosome * chromosome);

class Chromosome{
    protected:
        vector<Gene *>		*genes;
        double				fitness;
		double				accFitness;
		double				accNormalizedFitness;
        FitnessFunction		fitnessFunction;
		vector<Gene *>		*geneModel;
		

    public:
        Chromosome(vector<Gene *> * geneModel, FitnessFunction fitFunction);
		Chromosome(Chromosome * original);
        ~Chromosome();
        vector<Gene *> *getGenes();
		vector<Gene *> *getGeneModel();
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


Chromosome::Chromosome(vector<Gene *> * genModel, FitnessFunction fitFunction){		
	geneModel = genModel;	
	genes = new vector<Gene *>();
	
	for (unsigned int i = 0; i < genModel->size(); i++) {
		Gene * gModel = (*genModel)[i];
		Gene *newGene = new Gene(gModel);		
		switch (newGene->getDataType())
		{
			case (FLOAT): {
				uniform_real_distribution<float> randomF(newGene->getMinimumSeed().floatValue, newGene->getMaximumSeed().floatValue);
				newGene->setValueFloat(randomF(genGA));
				break;
			}
			case (DOUBLE): {
				uniform_real_distribution<double> randomD(newGene->getMinimumSeed().doubleValue, newGene->getMaximumSeed().doubleValue);
				newGene->setValueDouble(randomD(genGA));
				break;
			}
			case (INT8): {
				uniform_int_distribution<int16_t> randomI16(newGene->getMinimumSeed().int8Value, newGene->getMaximumSeed().int8Value);
				#pragma warning(push)
				#pragma warning(disable : 4244)
				/*
				Should never cause loss of data, since seed values are within the 8 bit range
				*/
				newGene->setValueUInt8(randomI16(genGA));
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
				uniform_int_distribution<uint16_t> randomI16(newGene->getMinimumSeed().uint8Value, newGene->getMaximumSeed().uint8Value);			
				#pragma warning(push)
				#pragma warning(disable : 4244)
				/*
				Should never cause loss of data, since seed values are within the 8 bit range
				*/

				newGene->setValueUInt8(randomI16(genGA));
				#pragma warning(pop)
				break;
			}
			case (INT16): {
				uniform_int_distribution<int16_t> randomI16(newGene->getMinimumSeed().int16Value, newGene->getMaximumSeed().int16Value);
				newGene->setValueInt16(randomI16(genGA));
				break;
			}
			case (UINT16): {
				uniform_int_distribution<uint16_t> randomI16(newGene->getMinimumSeed().uint16Value, newGene->getMaximumSeed().uint16Value);
				newGene->setValueUInt16(randomI16(genGA));
				break;
			}
			case (INT32): {
				uniform_int_distribution<int32_t> randomI32(newGene->getMinimumSeed().int32Value, newGene->getMaximumSeed().int32Value);
				newGene->setValueInt32(randomI32(genGA));
				break;
			}
			case (UINT32): {
				uniform_int_distribution<uint32_t> randomI32(newGene->getMinimumSeed().uint32Value, newGene->getMaximumSeed().uint32Value);
				newGene->setValueUInt32(randomI32(genGA));
				break;
			}
			case (INT64): {
				uniform_int_distribution<int64_t> randomI64(newGene->getMinimumSeed().int64Value, newGene->getMaximumSeed().int64Value);
				newGene->setValueInt64(randomI64(genGA));
				break;
			}
			case (UINT64): {
				uniform_int_distribution<uint64_t> randomI64(newGene->getMinimumSeed().uint64Value, newGene->getMaximumSeed().uint64Value);
				newGene->setValueUInt64(randomI64(genGA));
				break;
			}
			default: {
				//TODO: Set value for CUSTOM data type
				break;
			}
		}
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

	//Clone the vector
	genes = new vector<Gene *>(*original->getGenes());
}

Chromosome::~Chromosome(){
    delete genes;
}

vector<Gene *> * Chromosome::getGenes(){
    return genes;	
}

inline vector<Gene*> * Chromosome::getGeneModel()
{
	return geneModel;
}

inline void Chromosome::setChromosome(Chromosome &chm)
{
	vector<Gene *> * ind = chm.getGenes();
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
	for (vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {
		Gene * gen = *it;
		switch (gen->getDataType()) {
			case INT8:	
				printf("%6" PRIi8 , gen->getValueInt8());
				break;
			case UINT8:
				printf("%6" PRIu8, gen->getValueUInt8());
				break;
			case INT16:
				printf("%6" PRIi16, gen->getValueInt16());
				break;
			case UINT16:
				printf("%6" PRIu16, gen->getValueUInt16());
				break;
			case INT32:
				printf("%6" PRIi32, gen->getValueInt32());
			case UINT32:
				printf("%6" PRIu32, gen->getValueUInt32());
				break;
			case INT64:
				printf("%6" PRIi64, gen->getValueInt64());
			case UINT64:
				printf("%6" PRIu64, gen->getValueUInt64());
				break;
			case FLOAT:
				printf("%6.3f", gen->getValueFloat());
				break;
			case DOUBLE:			
				printf("%6.3f", gen->getValueDouble());
				break;
			case CUSTOM:
				std::cout << uppercase << hex << gen->getValueUInt64() << dec;
				break;			
		}		
		std::cout << ", \t";
	}
	std::cout << "F=" << fitness << endl;;
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
		if (
			((*genes)[i]->getDataType() == (*ind->getGenes())[i]->getDataType()) ||
			((*genes)[i]->getValueUInt64() != (*ind->getGenes())[i]->getValueUInt64())
			)
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

