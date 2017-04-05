#pragma once

#include "mutationvectorized.h"

class MutationUniform : public MutationVectorized {
	public:
		void		mutate(Chromosome *chromosome);
		void		mutate(Gene *gene);
};

inline void MutationUniform::mutate(Chromosome* chromosome)
{	
	std::vector<Gene *> *genes = chromosome->getGenes();
	for (std::vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {		
		mutate(*it);
	}
}

inline void MutationUniform::mutate(Gene* gene)
{
	gene->initialize();	
}