#pragma once

#define	PUNISHMENT_FACTOR		1e3f

#include "genetic-algorithm.h"

vector<Gene*> *getGenotype2d() {
	GeneValue minSeed, maxSeed, minSeed2, maxSeed2;
	vector<Gene *> *genotype = new vector<Gene *>();

	minSeed.floatValue = 0;
	maxSeed.floatValue = 1.4f;

	maxSeed2.uint8Value = 1;
	minSeed2.uint8Value = 0;

	Gene *gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	genotype->push_back(gene);

	gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype->push_back(gene);

	return genotype;
}

double fitnessFunction2d(Chromosome * chromosome) {
	//http://www.zweigmedia.com/RealWorld/simplex.html

	float x1 = (*chromosome->getGenes())[0]->getValue().floatValue;	
	uint8_t y = (*chromosome->getGenes())[1]->getValue().uint8Value;

	//We have to check if the float/double values are valid
	if (isnan(x1) || isinf(x1))
		return -INFINITY;

	double fitness = -(-y + 2.0f * x1 - log(x1/2.0f));
	double punish = 0;

	bool violated = false;
	float violations[3] = { 0, 0, 0 };

	/*
	For constraint evaluation, we check if their complement is true.
	If so, the constraint has been violated and we should punish the fitness proportionately.
	*/

	if (-x1 - log(x1 / 2.0f) + (float)y > 0) {
		violated = true;
		violations[0] = -x1 - log(x1 / 2.0f) + (float)y;
	}
	if (x1 > 1.4f) {
		violated = true;
		violations[1] = 1.4f - x1;
	}

	if (x1 < 0.5f) {
		violated = true;
		violations[2] = 0.5f - x1;
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <3; i++)
			punish -= abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}


vector<Gene*> *getGenotype2() {
	GeneValue minSeed, maxSeed, minSeed2, maxSeed2;
	vector<Gene *> *genotype = new vector<Gene *>();

	minSeed.floatValue = 0;
	maxSeed.floatValue = 1.4f;

	maxSeed2.uint8Value = 1;
	minSeed2.uint8Value = 0;

	Gene *gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	genotype->push_back(gene);

	gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	genotype->push_back(gene);

	gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype->push_back(gene);

	return genotype;
}


double fitnessFunction2(Chromosome * chromosome) {
	//http://www.zweigmedia.com/RealWorld/simplex.html

	float x1 = (*chromosome->getGenes())[0]->getValue().floatValue;
	float x2 = (*chromosome->getGenes())[1]->getValue().floatValue;
	uint8_t y = (*chromosome->getGenes())[2]->getValue().uint8Value;

	//We have to check if the float/double values are valid
	if (isnan(x1) || isinf(x1))
		return -INFINITY;
	if (isnan(x2) || isinf(x2))
		return -INFINITY;

	double fitness = -(-(float)y + 2.0f * x1 + x2);
	double punish = 0;

	bool violated = false;
	float violations[4] = { 0, 0, 0, 0 };

	/*
	For constraint evaluation, we check if their complement is true.
	If so, the constraint has been violated and we should punish the fitness proportionately.
	*/

	if (x1 - 2.0f * exp(-x2) != 0) {
		violated = true;
		violations[0] = x1 - 2.0f * exp(-x2);
	}
	if (-x1 + x2 + (float)y > 0) {
		violated = true;
		violations[1] = -x1 + x2 + (float)y;
	}
	if (x1 > 1.4f) {
		violated = true;
		violations[2] = 1.4f - x1;
	}

	if (x1 < 0.5f) {
		violated = true;
		violations[3] = 0.5f - x1;
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <4; i++)
			punish -= abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}

vector<Gene*> *getGenotype1() {
	GeneValue minSeed, maxSeed, minSeed2, maxSeed2;
	vector<Gene *> *genotype = new vector<Gene *>();

	minSeed.floatValue = 0;
	maxSeed.floatValue = 1.6f;

	maxSeed2.uint8Value = 1;
	minSeed2.uint8Value = 0;

	Gene *gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	genotype->push_back(gene);

	gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype->push_back(gene);

	return genotype;
}


double fitnessFunction1(Chromosome * chromosome) {
	//http://www.zweigmedia.com/RealWorld/simplex.html

	float x = (*chromosome->getGenes())[0]->getValue().floatValue;
	uint8_t y = (*chromosome->getGenes())[1]->getValue().uint8Value;

	//We have to check if the float/double values are valid
	if (isnan(x) || isinf(x))
		return -INFINITY;

	double fitness = -(2 * x + y);
	double punish = 0;

	bool violated = false;
	float violations[4] = { 0, 0, 0, 0 };

	if (1.25f > x*x + (float)y) {
		violated = true;
		violations[0] = 1.25f - (x*x + (float)y);
	}
	if (x + y > 1.6f) {
		violated = true;
		violations[1] = 1.6f - (x + y);
	}
	if (x > 1.6f) {
		violated = true;
		violations[2] = 1.6f - x;
	}

	if (x < 0) {
		violated = true;
		violations[3] = -x;
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <4; i++)
			punish -= abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}