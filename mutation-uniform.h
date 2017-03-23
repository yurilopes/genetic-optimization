#pragma once

#include "mutation.h"
#include "chromosome.h"

class MutationUniform : public Mutation {
	public:
		void		mutate(Chromosome *chromosome);
};

inline void MutationUniform::mutate(Chromosome* chromosome)
{	
	std::vector<Gene *> *genes = chromosome->getGenes();
	for (std::vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {
		Gene * gene = *it;
		gene->initialize();
	}
}