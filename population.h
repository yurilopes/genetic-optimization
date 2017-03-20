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
		vector<Gene *>						*genotype;
		void								deleteChromosomes();
		
    public:
        Population(vector<Gene *> *genModel);
        ~Population();
        void								printPopulation();
        void								calculateFitness();		
        void								initialize(uint32_t populationSize, FitnessFunction fitFunction);
		vector<Chromosome*>					*getChromosomes();
		Chromosome							*getEqualChromosome(Chromosome * ind);
		Chromosome							*getFittestChromosome();
		void								refreshFitnessFunction(FitnessFunction func);
		vector<Gene *>						*getGenotype();

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

inline vector<Gene*>* Population::getGenotype()
{
	return genotype;
}

Population::~Population(){
	deleteChromosomes();
}

inline void Population::deleteChromosomes()
{
	for (vector<Chromosome *>::iterator it = chromosomes.begin(); it != chromosomes.end(); it++) {
		Chromosome * chm = *it;
		delete chm;
	}	
}

void Population::initialize(uint32_t populationSize, FitnessFunction fitFunction){
	deleteChromosomes();
    for(uint32_t i=0; i<populationSize; i++){
        Chromosome *chm = new Chromosome(genotype, fitFunction);
        chromosomes.push_back(chm);
    }
}

Population::Population(vector<Gene *> *genModel){
	genotype = genModel;
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