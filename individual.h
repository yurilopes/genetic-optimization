#include <vector>

typedef int32_t (*FitnessFunction)(vector<int32_t>*);

class Individual{
    protected:
        vector<int32_t> *individual;
        int32_t         fitness;
        FitnessFunction fitnessFunction;

    public:
        Individual(int32_t indSize, FitnessFunction fitFunction);
        ~Individual();
        vector<int32_t>* getIndividual();
        void calculateFitness();
        int32_t getFitness();
};


Individual::Individual(int32_t indSize, FitnessFunction fitFunction){
    individual = new vector<int32_t>(indSize);
    fitness = 0;
    fitnessFunction = fitFunction;
}

Individual::~Individual(){
    delete individual;
}

vector<int32_t>* Individual::getIndividual(){
    return individual;
}

int32_t Individual::getFitness(){
    return fitness;
}

void Individual::calculateFitness(){
    fitness=fitnessFunction(individual);
}
