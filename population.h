#include <list>
#include <vector>
#include <random>
#include <ctime>

#include "individual.h"
#include "utils.h"

class Population{
    protected:
        list<Individual*> population;
        mt19937 *gen;
        uniform_int_distribution<int32_t> *disLong;
        int32_t random();
        void createPopulation(int32_t populationSize, int32_t individualSize, FitnessFunction fitFunction);

    public:
        Population(int32_t populationSize, int32_t individualSize, FitnessFunction fitFunction);
        ~Population();
        void printPopulation();
        void calculateFitness();

};

Population::~Population(){
    for(list<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
        delete individual;
    }
}

int32_t Population::random(){
    return (*disLong)((*gen));
}

void Population::createPopulation(int32_t populationSize, int32_t individualSize, FitnessFunction fitFunction){
    for(int32_t i=0; i<populationSize; i++){
        Individual *individual = new Individual(individualSize, fitFunction);
        for(int32_t j=0; j<individualSize; j++){
            (* individual->getIndividual())[j]=random();
        }
        population.push_back(individual);
    }
}

Population::Population(int32_t populationSize, int32_t individualSize, FitnessFunction fitFunction){

    gen = new mt19937(static_cast<unsigned int>(std::time(0)));
    disLong = new uniform_int_distribution<int32_t>(SHRT_MIN, SHRT_MAX); //(LONG_MIN, LONG_MAX);

    createPopulation(populationSize, individualSize, fitFunction);
}


void Population::printPopulation(){
    for(list<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
        for(vector<int32_t>::iterator it = individual->getIndividual()->begin(); it!=individual->getIndividual()->end(); it++){
            cout<<*it;
            cout<<", ";
        }
        cout<<"F="<<individual->getFitness();
        cout<<endl;
    }
}

void Population::calculateFitness(){
    for(list<Individual*>::iterator itr = population.begin(); itr!=population.end(); itr++){
        Individual *individual = *itr;
        individual->calculateFitness();
    }
}
