#pragma once

#include "mutation.h"

class MutationVectorized : public Mutation { //Abstract class
public:
	virtual void		mutate(Chromosome *chromosome) = 0; //pure virtual function
	virtual void		mutate(Gene *gene) = 0; //pure virtual function
};