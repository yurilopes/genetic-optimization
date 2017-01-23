#include <list>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <limits.h>

#include "individual.h"
#include "utils.h"

class Population{
    protected:
        list<Individual*> population;
        mt19937 *gen;
        uniform_int_distribution<int32_t> *dis;
        int32_t random();

    public:
        Population(int32_t minGeneValue, int32_t maxGeneValue);
        ~Population();
        void printPopulation();
        void calculateFitness();
		void calculateSelection();
        void initialize(int32_t populationSize, int32_t individualSize, FitnessFunction fitFunction);

};

Population::~Population(){
    for(list<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
        delete individual;
    }
}

int32_t Population::random(){
    return (*dis)((*gen));
}

void Population::initialize(int32_t populationSize, int32_t individualSize, FitnessFunction fitFunction){
    for(int32_t i=0; i<populationSize; i++){
        Individual *individual = new Individual(individualSize, fitFunction);
        for(int32_t j=0; j<individualSize; j++){
            (* individual->getIndividual())[j]=random();
        }
        population.push_back(individual);
    }
}

Population::Population(int32_t minGeneValue, int32_t maxGeneValue){

    gen = new mt19937(static_cast<unsigned int>(std::time(0)));
    dis = new uniform_int_distribution<int32_t>(minGeneValue, maxGeneValue);

}


void Population::printPopulation(){
    for(list<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
        for(vector<int32_t>::iterator it = individual->getIndividual()->begin(); it!=individual->getIndividual()->end(); it++){
            cout<<*it;
            cout<<", ";
        }
        cout << "F=" << individual->getFitness();
		cout << ", FN=" << individual->getNormalizedFitness();
        cout<<endl;
    }
}


void Population::calculateFitness(){	
	/*
	Calculate fitness values for each individual
	*/
    for(list<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
        individual->calculateFitness();
    }

    population.sort( Individual::compare );
}


void Population::calculateSelection() {
	int32_t total = 0, lowest, current;
	/*
	Calculate the lowest fitness value
	This is used to translate every fitness value to a positive range
	*/
	for (list<Individual*>::iterator itr = population.begin(); itr != population.end(); itr++) {
		Individual *individual = *itr;
		current = individual->getFitness();
		if (itr == population.begin()) //Lowest = fitness of the first iteration
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
	for (list<Individual*>::iterator itr = population.begin(); itr != population.end(); itr++) {
		Individual *individual = *itr;
		total+=individual->addToFitness(lowest);
	}

	/*
	Calculate the normalized fitness for each individual
	*/
	for (list<Individual*>::iterator itr = population.begin(); itr != population.end(); itr++) {
		Individual *individual = *itr;
		individual->setNormalizedFitness((float)individual->getFitness() / (float)total);
	}

	

}
