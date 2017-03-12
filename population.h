#pragma once

#include <vector>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "chromosome.h"
#include "utils.h"

class Population{
    protected:
        vector<Chromosome*>					chromosomes;
        mt19937								*gen;
        uniform_int_distribution<int32_t>	*dis;
        int32_t								random();		
		vector<Gene *>						*geneModel;
		

    public:
        Population(vector<Gene *> *genModel);
        ~Population();
        void printPopulation();
        void calculateFitness();		
        void initialize(uint32_t populationSize, FitnessFunction fitFunction);
		vector<Chromosome*>* getChromosomes();
		Chromosome * getEqualChromosome(Chromosome * ind);
		Chromosome * getFittestChromosome();
		void refreshFitnessFunction(FitnessFunction func);

};

vector<Chromosome*>* Population::getChromosomes() {
	return &chromosomes;
}

inline Chromosome * Population::getEqualChromosome(Chromosome * ind)
{
	for (vector<Chromosome*>::iterator itr = chromosomes.begin(); itr != chromosomes.end(); itr++) {
		Chromosome *individual = *itr;
		if (individual->equals(ind))
			return individual;
	}
	return NULL;
}

inline Chromosome * Population::getFittestChromosome()
{
	return chromosomes.front();
}

inline void Population::refreshFitnessFunction(FitnessFunction func)
{
	for (vector<Chromosome*>::iterator itr = chromosomes.begin(); itr != chromosomes.end(); itr++) {
		Chromosome *individual = *itr;
		individual->setFitnessFunction(func);
	}
}

Population::~Population(){
	for(vector<Chromosome*>::iterator itr = chromosomes.begin(); itr!=chromosomes.end(); itr++){
		Chromosome *individual = *itr;
		delete individual;
	}
}

int32_t Population::random(){
    return (*dis)((*gen));
}

void Population::initialize(uint32_t populationSize, FitnessFunction fitFunction){
    for(uint32_t i=0; i<populationSize; i++){
        Chromosome *individual = new Chromosome(geneModel, fitFunction);
        chromosomes.push_back(individual);
    }
}

Population::Population(vector<Gene *> *genModel){
	geneModel = genModel;
}


void Population::printPopulation(){
	uint32_t i = 0;
    for(vector<Chromosome*>::iterator itr = chromosomes.begin(); itr!=chromosomes.end(); itr++){
        Chromosome *chrm = *itr;
		cout << i << ": ";
		chrm->print();
		//cout << "\t, FN=" << individual->getAccNormalizedFitness();
        //cout<<endl;
		i++;
    }
}


void Population::calculateFitness(){	
	/*
	Calculate fitness values for each chromosome
	*/
    for(vector<Chromosome*>::iterator itr = chromosomes.begin(); itr!=chromosomes.end(); itr++){
        Chromosome *chm = *itr;
        chm->calculateFitness();
    }

    sort(chromosomes.begin(), chromosomes.end(), Chromosome::compare);
}