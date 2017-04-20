#pragma once

#include "mutationvectorized.h"

class MutationVector : public Mutation {
	protected:
		std::vector<MutationVectorized *> *mutationVector = NULL;

	public:
		void		mutate(Chromosome *chromosome);
		void		setMutationVector(std::vector<MutationVectorized *> mutVector);
		std::vector<MutationVectorized *> *getMutationVector();

		MutationVector(std::vector<MutationVectorized *> mutVector);
		~MutationVector();
};

inline void MutationVector::mutate(Chromosome *chromosome) {
	if (!mutationVector)
		return;
	
	size_t i = 0;

	std::vector<Gene *> *genes = chromosome->getGenes();
	for (std::vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {
		if (i >= mutationVector->size())
			return;
		if((*mutationVector)[i]) //Ignore null pointers
			(*mutationVector)[i]->mutate(*it);				

		i++;
	}
}

inline void MutationVector::setMutationVector(std::vector<MutationVectorized*> mutVector)
{
	if (mutationVector)
		delete mutationVector;

	mutationVector = new std::vector<MutationVectorized*>();
	for (std::vector<MutationVectorized *>::iterator it = mutVector.begin(); it != mutVector.end(); it++) {
		mutationVector->push_back(*it);
	}
}

inline std::vector<MutationVectorized*>* MutationVector::getMutationVector()
{
	return mutationVector;
}

inline MutationVector::MutationVector(std::vector<MutationVectorized*> mutVector)
{
	mutationVector = new std::vector<MutationVectorized*>();
	for (std::vector<MutationVectorized *>::iterator it = mutVector.begin(); it != mutVector.end(); it++) {
		mutationVector->push_back(*it);
	}
}

inline MutationVector::~MutationVector()
{
	if (mutationVector)
		delete mutationVector;
}
