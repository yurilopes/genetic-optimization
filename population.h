#include <vector>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "individual.h"
#include "utils.h"

class Population{
    protected:
        vector<Individual*> population;
        mt19937 *gen;
        uniform_int_distribution<int32_t> *dis;
        int32_t random();
		bool	mainPopulation = true;
		

    public:
        Population(int32_t minSeed, int32_t maxSeed);
        ~Population();
        void printPopulation();
        void calculateFitness();		
        void initialize(uint32_t populationSize, uint32_t individualSize, FitnessFunction fitFunction);
		vector<Individual*>* getIndividualVector();
		Individual * getFittestIndividual();
		void refreshFitnessFunction(FitnessFunction func);
		bool isMainPopulation();
		void setMainPopulation(bool main);

};

vector<Individual*>* Population::getIndividualVector() {
	return &population;
}

inline Individual * Population::getFittestIndividual()
{
	return population.front();
}

inline void Population::refreshFitnessFunction(FitnessFunction func)
{
	for (vector<Individual*>::iterator itr = population.begin(); itr != population.end(); itr++) {
		Individual *individual = *itr;
		individual->setFitnessFunction(func);
	}
}

inline bool Population::isMainPopulation()
{
	return mainPopulation;
}

inline void Population::setMainPopulation(bool main)
{
	mainPopulation = main;
}

Population::~Population(){
	if(mainPopulation) //The destructor won't try to free invalid Individual pointers
		for(vector<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
			Individual *individual = *itr;
			delete individual;
		}
}

int32_t Population::random(){
    return (*dis)((*gen));
}

void Population::initialize(uint32_t populationSize, uint32_t individualSize, FitnessFunction fitFunction){
    for(uint32_t i=0; i<populationSize; i++){
        Individual *individual = new Individual(individualSize, fitFunction);
        for(uint32_t j=0; j<individualSize; j++){
            (* individual->getIndividual())[j]=random();
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
    for(vector<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
		cout << i << ": ";
        for(vector<int32_t>::iterator it = individual->getIndividual()->begin(); it!=individual->getIndividual()->end(); it++){
            cout << setw(6) << *it;
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
    for(vector<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
        individual->calculateFitness();
    }

    sort(population.begin(), population.end(), Individual::compare);
}