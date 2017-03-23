#pragma once

#include "chromosome.h"

class Crossover { //Abstract class
public:
	virtual void		crossover(Chromosome * parentX, Chromosome * parentY, Chromosome * childX, Chromosome * childY) = 0; //pure virtual function
};