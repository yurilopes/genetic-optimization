#pragma once

#include "crossover.h"
#include "genetic-algorithm.h"
#include <bitset>

class CrossoverUniformBitwise : public Crossover {
	void		crossover(Chromosome * parentX, Chromosome * parentY, Chromosome * childX, Chromosome * childY);
};

void CrossoverUniformBitwise::crossover(Chromosome * parentX, Chromosome * parentY, Chromosome * childX, Chromosome * childY) {
	/*
	Uniform crossover
	More info at: https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)

	This function changes both parents to the result of the crossover operation
	*/

	/*
									****************
									*** WARNING! ***
									****************

	Treating the GeneValue structure as a uint64_t should be fine
	double and uint64_t should have the same size (64 bits, 8 bytes) and they're the biggest data types in GeneValue
	*/

	bool first = false;

	for (size_t i = 0; i < parentX->getGenes()->size(); i++) {
		//First we iterate through the genes

		Gene *genXV = (*parentX->getGenes())[i];
		Gene *genYV = (*parentY->getGenes())[i];		

		bitset<sizeof(uint64_t)*CHAR_BIT> bX(genXV->getValue().uint64Value);
		bitset<sizeof(uint64_t)*CHAR_BIT> bY(genYV->getValue().uint64Value);

		for (size_t j = 0; j < sizeof(uint64_t)*CHAR_BIT; j++) {
			//Then we iterate through the bits of each gene pair
			if (i == parentX->getGenes()->size() - 1 && j == sizeof(uint64_t)*CHAR_BIT - 1) {
				/*
				The first bit is always kept untouched
				This is a premise of crossover
				*/
				continue;
			}

			if (randomBinary()) {
				/*
				Set bitX[j] as bitY[j] and vice-versa
				*/
				bool bitX = bX.test(j);
				bX[j] = bY[j];
				bY[j] = bitX;
			}
		}

		//Replace the children's old values with the new ones
		GeneValue value;
		value.uint64Value = bX.to_ullong();
		(*childX->getGenes())[i]->setValue(value);
		value.uint64Value = bY.to_ullong();
		(*childY->getGenes())[i]->setValue(value);
	}
}