//Comment this line below if running in Visual Studio
#include "stdafx.h"

#include <limits.h>

using namespace std;

#include "genetic-algorithm.h"
#include "population.h"


int32_t fitnessFunction(vector<int32_t> *variables){
    int32_t soma = 0;
    for(vector<int32_t>::iterator it = variables->begin(); it!=variables->end(); it++){
        soma += *it;
    }
    return soma;
}

int main(){
    /*
    A population is list of individuals
    Each individual is vector of ints
    Each problem variable/gene is an int
    */
    //Population population(LONG_MIN, LONG_MAX);
	Population population(-5, 5);
    population.initialize(5, 3, fitnessFunction);
    population.calculateFitness();
	population.calculateSelection();
	population.printPopulation();

    return 0;
}
