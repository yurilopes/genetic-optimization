/*
Uniform crossover
More info at: https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)

*/

#pragma once

#include "crossover.h"
#include <random>
#include <bitset>
#include <climits>

class CrossoverUniformBitwise : public Crossover {
	protected:
		std::uniform_int_distribution<uint16_t>		*distribution = NULL;
		std::mt19937								*randomGen = NULL;

	public:
		void		crossover(Chromosome * parentX, Chromosome * parentY, Chromosome * childX, Chromosome * childY);

		CrossoverUniformBitwise();
		~CrossoverUniformBitwise();
};

void CrossoverUniformBitwise::crossover(Chromosome * parentX, Chromosome * parentY, Chromosome * childX, Chromosome * childY) {

	/*
									****************
									*** WARNING! ***
									****************

	Treating the GeneValue structure as a uint64_t should be fine
	double and uint64_t should have the same size (64 bits, 8 bytes) and they're the biggest data types in GeneValue
	*/

	for (size_t i = 0; i < parentX->getGenes()->size(); i++) {
		//First we iterate through the genes

		Gene *genXV = (*parentX->getGenes())[i];
		Gene *genYV = (*parentY->getGenes())[i];

		bitset<sizeof(uint64_t)*CHAR_BIT> bX(genXV->getValue().uint64Value);
		bitset<sizeof(uint64_t)*CHAR_BIT> bY(genYV->getValue().uint64Value);

		for (size_t j = 0; j < sizeof(uint64_t)*CHAR_BIT; j++) {
			//Then we iterate through the bits of each gene pair
		//	if (i == 0 && j == sizeof(uint64_t)*CHAR_BIT - 1) {
				/*
				The first bit is always kept untouched
				This is a premise of crossover
				*/
		//		continue;
		//	}

			if ((*distribution)(*randomGen)) { //Generates 0 or 1 uniformly distributed
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

inline CrossoverUniformBitwise::CrossoverUniformBitwise()
{
	randomGen = new std::mt19937(static_cast<unsigned int>(std::time(0)));
	distribution = new std::uniform_int_distribution<uint16_t>(0, 1);
}

inline CrossoverUniformBitwise::~CrossoverUniformBitwise()
{
	if (distribution)
		delete distribution;
	if (randomGen)
		delete randomGen;
}
