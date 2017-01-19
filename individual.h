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

        static bool compare(Individual *ind0, Individual *ind1);
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

bool Individual::compare(Individual *ind0, Individual *ind1){
    /*
    Minimization of the fitness function is implied here
    */
    return ind0->getFitness() < ind1->getFitness();
}

