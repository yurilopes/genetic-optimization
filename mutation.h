#pragma once

#include "chromosome.h"

class Mutation { //Abstract class
	public:
		virtual void		mutate(Chromosome *chromosome) = 0; //pure virtual function
};