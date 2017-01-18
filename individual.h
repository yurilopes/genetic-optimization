#include <vector>

class Individual{
    protected:
        vector<int32_t> *individual;
        int32_t         fitness;

    public:
        Individual(int32_t indSize);
        ~Individual();
        vector<int32_t>* getIndividual();
        int32_t getFitness();
};


Individual::Individual(int32_t indSize){
    individual = new vector<int32_t>(indSize);
    fitness = 0;
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
