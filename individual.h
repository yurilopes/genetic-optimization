#include <vector>
#include <conio.h>

typedef int32_t (*FitnessFunction)(vector<int32_t>*);

class Individual{
    protected:
        vector<int32_t> *individual;
        int32_t         fitness;
		uint32_t        accFitness;
		float			accNormalizedFitness;
        FitnessFunction fitnessFunction;

    public:
        Individual(int32_t indSize, FitnessFunction fitFunction);
		Individual(Individual * original);
        ~Individual();
        vector<int32_t>* getIndividual();
		void setIndividual(vector<int32_t>* ind);
		int32_t calculateFitness();
		int32_t addToAccFitness(int32_t value);
        int32_t getFitness();
		uint32_t getAccFitness();
		void setFitnessFunction(FitnessFunction fitFunc);
		FitnessFunction getFitnessFunction();
		float getAccNormalizedFitness();
		void setAccNormalizedFitness(float fit);
		void print();

        static bool compare(Individual *ind0, Individual *ind1);
};


Individual::Individual(int32_t indSize, FitnessFunction fitFunction){
    individual = new vector<int32_t>(indSize);
    fitness = 0;
	accNormalizedFitness = 0.0;
    fitnessFunction = fitFunction;
}

inline Individual::Individual(Individual * original)
{
	fitness = original->getFitness();
	accNormalizedFitness = original->getAccNormalizedFitness();
	fitnessFunction = original->getFitnessFunction();

	//Clone the vector
	individual = new vector<int32_t>(*original->getIndividual());
}

Individual::~Individual(){
    delete individual;
}

vector<int32_t>* Individual::getIndividual(){
    return individual;
}

inline void Individual::setIndividual(vector<int32_t>* ind)
{
	if (individual != NULL)
		delete individual;
	individual = ind;
}

int32_t Individual::getFitness(){
    return fitness;
}

inline uint32_t Individual::getAccFitness()
{
	return accFitness;
}

inline void Individual::setFitnessFunction(FitnessFunction fitFunc)
{
	fitnessFunction = fitFunc;
}

inline FitnessFunction Individual::getFitnessFunction()
{
	return fitnessFunction;
}

float Individual::getAccNormalizedFitness() {
	return accNormalizedFitness;
}

void Individual::setAccNormalizedFitness(float fit)
{
	accNormalizedFitness = fit;
}

inline void Individual::print()
{
	for (vector<int32_t>::iterator it = individual->begin(); it != individual->end(); it++) {
		cout << left << setw(6) << *it;
		cout << ", \t";
	}
	cout << "F=" << fitness;	
	cout << endl;
}

int32_t Individual::calculateFitness(){
	fitness = fitnessFunction(individual);
	accFitness = fitness;
    return fitness;
}

inline int32_t Individual::addToAccFitness(int32_t value)
{
	return accFitness+=value;
}

bool Individual::compare(Individual *ind0, Individual *ind1){
    /*
    Maximization of the fitness function is implied here
    */
    return ind0->getFitness() > ind1->getFitness();
}

