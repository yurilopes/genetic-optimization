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
        vector<Chromosome*> population;
        mt19937 *gen;
        uniform_int_distribution<int32_t> *dis;
        int32_t random();		
		

    public:
        Population(int32_t minSeed, int32_t maxSeed);
        ~Population();
        void printPopulation();
        void calculateFitness();		
        void initialize(uint32_t populationSize, uint32_t individualSize, FitnessFunction fitFunction);
		vector<Chromosome*>* getIndividualVector();
		Chromosome * getEqualIndividual(Chromosome * ind);
		Chromosome * getFittestIndividual();
		void refreshFitnessFunction(FitnessFunction func);

};

vector<Chromosome*>* Population::getIndividualVector() {
	return &population;
}

inline Chromosome * Population::getEqualIndividual(Chromosome * ind)
{
	for (vector<Chromosome*>::iterator itr = population.begin(); itr != population.end(); itr++) {
		Chromosome *individual = *itr;
		if (individual->equals(ind))
			return individual;
	}
	return NULL;
}

inline Chromosome * Population::getFittestIndividual()
{
	return population.front();
}

inline void Population::refreshFitnessFunction(FitnessFunction func)
{
	for (vector<Chromosome*>::iterator itr = population.begin(); itr != population.end(); itr++) {
		Chromosome *individual = *itr;
		individual->setFitnessFunction(func);
	}
}

Population::~Population(){
	for(vector<Chromosome*>::iterator itr = population.begin(); itr!=population.end(); itr++){
		Chromosome *individual = *itr;
		delete individual;
	}
}

int32_t Population::random(){
    return (*dis)((*gen));
}

void Population::initialize(uint32_t populationSize, uint32_t individualSize, FitnessFunction fitFunction){
    for(uint32_t i=0; i<populationSize; i++){
        Chromosome *individual = new Chromosome(individualSize, fitFunction);
        for(uint32_t j=0; j<individualSize; j++){
            (* individual->getGenes())[j]=random();
        }
        population.push_back(individual);
    }
}

Population::Population(int32_t minSeed, int32_t maxSeed){

    gen = new mt19937(static_cast<unsigned int>(std::time(0)));
    dis = new uniform_int_distribution<int32_t>(minSeed, maxSeed);

}


void Population::printPopulation(){
	uint32_t i = 0;
    for(vector<Chromosome*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Chromosome *individual = *itr;
		cout << i << ": ";
        for(vector<int32_t>::iterator it = individual->getGenes()->begin(); it!=individual->getGenes()->end(); it++){
            cout << left << setw(6) << *it;
            cout << ", \t";
        }
        cout << "F=" << individual->getFitness();
		cout << "\t, FN=" << individual->getAccNormalizedFitness();
        cout<<endl;
		i++;
    }
}


void Population::calculateFitness(){	
	/*
	Calculate fitness values for each individual
	*/
    for(vector<Chromosome*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Chromosome *individual = *itr;
        individual->calculateFitness();
    }

    sort(population.begin(), population.end(), Chromosome::compare);
}