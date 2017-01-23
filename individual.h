#include <vector>

typedef int32_t (*FitnessFunction)(vector<int32_t>*);

class Individual{
    protected:
        vector<int32_t> *individual;
        int32_t         fitness;
		float			normalizedFitness;
        FitnessFunction fitnessFunction;

    public:
        Individual(int32_t indSize, FitnessFunction fitFunction);
        ~Individual();
        vector<int32_t>* getIndividual();
		int32_t calculateFitness();
		int32_t addToFitness(int32_t value);
        int32_t getFitness();
		float getNormalizedFitness();
		void setNormalizedFitness(float fit);

        static bool compare(Individual *ind0, Individual *ind1);
};


Individual::Individual(int32_t indSize, FitnessFunction fitFunction){
    individual = new vector<int32_t>(indSize);
    fitness = 0;
	normalizedFitness = 0.0;
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

float Individual::getNormalizedFitness() {
	return normalizedFitness;
}

void Individual::setNormalizedFitness(float fit)
{
	normalizedFitness = fit;
}

int32_t Individual::calculateFitness(){
    return fitness=fitnessFunction(individual);
}

inline int32_t Individual::addToFitness(int32_t value)
{
	return fitness+=value;
}

bool Individual::compare(Individual *ind0, Individual *ind1){
    /*
    Maximization of the fitness function is implied here
    */
    return ind0->getFitness() > ind1->getFitness();
}

