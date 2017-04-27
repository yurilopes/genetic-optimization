#pragma once

#define	PUNISHMENT_FACTOR		-1e3f

#include "genetic-algorithm.h"

int its = 0;

vector<Gene*> *getGenotype7() {
	GeneValue minSeed, maxSeed;
	vector<Gene *> *genotype = new vector<Gene *>();

	Gene * gene;

	minSeed.uint8Value = 1;
	maxSeed.uint8Value = 3;

	//N1..N3
	for (int i = 0; i < 3; i++) {
		gene = new Gene(UINT8);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	minSeed.floatValue = 250.0f;
	maxSeed.floatValue = 2500.0f;

	//V1..V3
	for (int i = 0; i < 3; i++) {		
		gene = new Gene(FLOAT);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	//B1..B2
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			//B1
			minSeed.floatValue = 400.0f / 9.0f; //44.4f
			maxSeed.floatValue = 625.0f;
		}
		if (i == 1) {
			//B2
			minSeed.floatValue = 160.0f / 9.0f; //17.7f
			maxSeed.floatValue = 1250.0f / 3.0f; //416.6f
		}
		gene = new Gene(FLOAT);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	//TL1..TL2
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			//TL1
			minSeed.floatValue = 20.0f / 3.0f; //6.66f
			maxSeed.floatValue = 20.0f;
		}
		if (i == 1) {
			//TL2
			minSeed.floatValue = 16.0f / 3.0f; //5.33f
			maxSeed.floatValue = 16.0f;
		}
		gene = new Gene(FLOAT);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	return genotype;
}

double fitnessFunction7(Chromosome * chromosome) {

	double N[3];
	double V[3];
	double B[2];
	double TL[2];

	N[0] = (double)(*chromosome->getGenes())[0]->getValue().uint8Value;
	N[1] = (double)(*chromosome->getGenes())[1]->getValue().uint8Value;
	N[2] = (double)(*chromosome->getGenes())[2]->getValue().uint8Value;
	
	V[0] = (double)(*chromosome->getGenes())[3]->getValue().floatValue;
	V[1] = (double)(*chromosome->getGenes())[4]->getValue().floatValue;
	V[2] = (double)(*chromosome->getGenes())[5]->getValue().floatValue;

	B[0] = (double)(*chromosome->getGenes())[6]->getValue().floatValue;
	B[1] = (double)(*chromosome->getGenes())[7]->getValue().floatValue;

	TL[0] = (double)(*chromosome->getGenes())[8]->getValue().floatValue;
	TL[1] = (double)(*chromosome->getGenes())[9]->getValue().floatValue;

	double S[2][3] = {
		{ 2.0, 3.0, 4.0 },
		{ 4.0, 6.0, 3.0 }
	};

	double t[2][3] = {
		{ 8.0, 20.0, 8.0 },
		{ 16.0, 4.0, 4.0 }
	};

	double Q[2] = { 40000.0, 20000.0 };

	double fitness = 0;
	for (int j = 0; j < 3; j++)
		fitness -= 250.0 * N[j] * pow(V[j], 0.6);

	double punish = 0;

	bool violated = false;
	double violations[14] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	int iv = 0;

	for (int i = 0; i < 2; i++) {
		if ((Q[i] * TL[i]) / B[i] > 6000.0) {
			violated = true;
			violations[iv] = 6000.0 - ((Q[i] * TL[i]) / B[i]);
		}
		iv++;
	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			if (V[j] < S[i][j] * B[i]) {
				violated = true;
				violations[iv] = S[i][j] * B[i] - (V[j]);
			}
			iv++;
		}
	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			if (N[j] * TL[i] < t[i][j]) {
				violated = true;
				violations[iv] = t[i][j] - (N[j] * TL[i]);
			}
			iv++;
		}
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <14; i++)
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	if (its == 50) {
		cout << "--------------------------------------" << endl;
		chromosome->print();

		printf("iv=%d\nF=%f\nP=%f\n", iv, fitness, fitness + punish);
		system("pause");
	}

	return fitness + punish;
}

vector<Gene*> *getGenotype6() {
	GeneValue minSeed, maxSeed, minSeed2, maxSeed2;
	vector<Gene *> *genotype = new vector<Gene *>();

	minSeed.floatValue = 27.0f;
	maxSeed.floatValue = 45.0f;

	Gene * gene;

	//x1..x3
	for (int i = 0; i < 3; i++) {
		gene = new Gene(FLOAT);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	//y1	
	minSeed2.uint8Value = 78;
	maxSeed2.uint8Value = 102;
	gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype->push_back(gene);

	//y2
	minSeed2.uint8Value = 33;
	maxSeed2.uint8Value = 45;
	gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype->push_back(gene);

	return genotype;
}

double fitnessFunction6(Chromosome * chromosome) {

	double x1 = (double)(*chromosome->getGenes())[0]->getValue().floatValue;
	double x2 = (double)(*chromosome->getGenes())[1]->getValue().floatValue;
	double x3 = (double)(*chromosome->getGenes())[2]->getValue().floatValue;
	double y1 = (double)(*chromosome->getGenes())[3]->getValue().uint8Value;
	double y2 = (double)(*chromosome->getGenes())[4]->getValue().uint8Value;

	//We have to check if the float/double values are valid
	if (isnan(x1) || isinf(x1))
		return -INFINITY;
	if (isnan(x2) || isinf(x2))
		return -INFINITY;
	if (isnan(x3) || isinf(x3))
		return -INFINITY;

	double fitness = -5.357854*x1*x1 - 0.835689*y1*x3 - 37.29329*y1 + 40792.141;
	double punish = 0;

	bool violated = false;
	double violations[3] = { 0, 0, 0 };

	double s1 = 85.334407 + 0.0056858*y2*x3 + 0.0006262*y1*x2 - 0.0022053*x1*x3;
	double s2 = 80.51249 + 0.0071317*y2*x3 + 0.0029955*y1*y2 + 0.0021813*pow(x1, 2.0) - 90.0;
	double s3 = 9.300961 + 0.0047026*x1*x3 + 0.0012547*y1*x1 + 0.0019085*x1*x2 - 20.0;

	if (s1 > 92.0) {
		violated = true;
		violations[0] = 92.0 - s1;		
	}
	if (s2 > 20.0) {
		violated = true;
		violations[1] = 20.0 - s2;		
	}
	if (s3 > 5.0) {
		violated = true;
		violations[2] = 5.0 - s3;		
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <3; i++)
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}

vector<Gene*> *getGenotype5() {
	GeneValue minSeed, maxSeed, minSeed2, maxSeed2;
	vector<Gene *> *genotype = new vector<Gene *>();

	minSeed.doubleValue = 0.0f;
	maxSeed.doubleValue = 5.0f;

	minSeed2.uint8Value = 0;
	maxSeed2.uint8Value = 1;
	
	Gene * gene;

	//x1..x3
	for (int i = 0; i < 3; i++) {
		gene = new Gene(DOUBLE);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}
		
	//y1..y4
	for (int i = 0; i < 4; i++) {
		gene = new Gene(UINT8);
		gene->setSeedRange(minSeed2, maxSeed2);
		gene->enableBounding(true);
		gene->setBounds(minSeed2, maxSeed2);
		genotype->push_back(gene);
	}

	return genotype;
}

double fitnessFunction5(Chromosome * chromosome) {	

	double x1 = (*chromosome->getGenes())[0]->getValue().doubleValue;
	double x2 = (*chromosome->getGenes())[1]->getValue().doubleValue;
	double x3 = (*chromosome->getGenes())[2]->getValue().doubleValue;
	double y1 = (double)(*chromosome->getGenes())[3]->getValue().uint8Value;
	double y2 = (double)(*chromosome->getGenes())[4]->getValue().uint8Value;
	double y3 = (double)(*chromosome->getGenes())[5]->getValue().uint8Value;
	double y4 = (double)(*chromosome->getGenes())[6]->getValue().uint8Value;

	//We have to check if the float/double values are valid
	if (isnan(x1) || isinf(x1))
		return -INFINITY;
	if (isnan(x2) || isinf(x2))
		return -INFINITY;
	if (isnan(x3) || isinf(x3))
		return -INFINITY;

	double fitness = -(pow((y1 - 1.0), 2.0) + pow((y2 - 2.0), 2.0) + pow((y3 - 1.0), 2.0) - log(y4 + 1.0) + pow((x1 - 1.0), 2.0) + pow((x2 - 2.0), 2.0) + pow((x3 - 3.0), 2.0));
	double punish = 0;

	bool violated = false;
	double violations[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	if (y1 + y2 + y3 + x1 + x2 + x3 > 5.0) {
		violated = true;
		violations[0] = 5.0 - (y1 + y2 + y3 + x1 + x2 + x3);
	}
	if (pow(y3, 2.0) + pow(x1, 2.0) + pow(x2, 2.0) + pow(x3, 2.0) > 5.5) {
		violated = true;
		violations[1] = 5.5 - (pow(y3, 2.0) + pow(x1, 2.0) + pow(x2, 2.0) + pow(x3, 2.0));
	}
	if (y1 + x1 > 1.2) {
		violated = true;
		violations[2] = 1.2 - (y1 + x1);
	}
	if (y2 + x2 > 1.8) {
		violated = true;
		violations[3] = 1.8 - (y2 + x2);
	}
	if (y3 + x3 > 2.5) {
		violated = true;
		violations[4] = 2.5 - (y3 + x3);
	}
	if (y4 + x1 > 1.2) {
		violated = true;
		violations[5] = 1.2 - (y4 + x1);
	}
	if (pow(y2, 2.0) + pow(x2, 2.0) > 1.64) {
		violated = true;
		violations[6] = 1.64 - (pow(y2, 2.0) + pow(x2, 2.0));
	}
	if (pow(y3, 2.0) + pow(x3, 2.0) > 4.25) {
		violated = true;
		violations[7] = 4.25 - (pow(y3, 2.0) + pow(x3, 2.0));
	}
	if (pow(y2, 2.0) + pow(x3, 2.0) > 4.64) {
		violated = true;
		violations[8] = 4.64 - (pow(y2, 2.0) + pow(x3, 2.0));
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <9; i++)
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}

vector<Gene*> *getGenotype4() {
	GeneValue minSeed, maxSeed, minSeed2, maxSeed2;
	vector<Gene *> *genotype = new vector<Gene *>();	

	minSeed2.uint8Value = 0;
	maxSeed2.uint8Value = 1;

	//y
	Gene *gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype->push_back(gene);

	minSeed.floatValue = 0;
	maxSeed.floatValue = 10.0f; //3.514237f;

	//v1
	gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	gene->enableBounding(true);
	gene->setBounds(minSeed, maxSeed);
	genotype->push_back(gene);

	minSeed.floatValue = 0;
	maxSeed.floatValue = 10.0f;

	//v2
	gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	gene->enableBounding(true);
	gene->setBounds(minSeed, maxSeed);
	genotype->push_back(gene);

	minSeed.floatValue = 0;
	maxSeed.floatValue = 20.0f;

	//x1
	gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	gene->enableBounding(true);
	gene->setBounds(minSeed, maxSeed);
	genotype->push_back(gene);

	minSeed.floatValue = 0;
	maxSeed.floatValue = 20.0f;

	//x2
	gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	gene->enableBounding(true);
	gene->setBounds(minSeed, maxSeed);
	genotype->push_back(gene);

	return genotype;
}

double fitnessFunction4(Chromosome * chromosome) {
	//MODE_MAXIMIZATION

	float y = (float)(*chromosome->getGenes())[0]->getValue().uint8Value;
	float v1 = (*chromosome->getGenes())[1]->getValue().floatValue;
	float v2 = (*chromosome->getGenes())[2]->getValue().floatValue;
	float x1 = (*chromosome->getGenes())[3]->getValue().floatValue;
	float x2 = (*chromosome->getGenes())[4]->getValue().floatValue;

	//We have to check if the float/double values are valid
	if (isnan(v1) || isinf(v1))
		return INFINITY;
	if (isnan(v2) || isinf(v2))
		return INFINITY;
	if (isnan(x1) || isinf(x1))
		return INFINITY;
	if (isnan(x2) || isinf(x2))
		return INFINITY;

	float y2 = (1.0f - y);
	float zs = (0.9f*(1.0f - exp(-0.5f*v1)))*x1 + (0.8f*(1.0f - exp(-0.4f*v2)))*x2;

	double fitness = -(7.5f*y + 5.5f*y2 + 7.0f*v1 + 6.0f*v2 + 5.0f*(x1 + x2) );
	double punish = 0;	

	bool violated = false;
	float violations[5] = { 0, 0, 0, 0, 0 };

	if (zs != 10) {
		violated = true;
		violations[0] = 10 - zs;
	}
	if (v1 > 10*y) {
		violated = true;
		violations[1] = 10*y - v1;
	}
	if (v2 > 10 * y2) {
		violated = true;
		violations[2] = 10.0f*y2 - v2;		
	}
	if (x1 > 20.0f*y) {
		violated = true;
		violations[3] = 20.0f*y - x1;	
	}
	if (x2 > 20.0f*y2) {
		violated = true;
		violations[4] = 20.0f*y2 - x2;
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <5; i++)
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}

vector<Gene*> *getGenotype3() {
	//MODE_MAXIMIZATION
	GeneValue minSeed, maxSeed, minSeed2, maxSeed2, minSeed3, maxSeed3;
	vector<Gene *> *genotype = new vector<Gene *>();

	minSeed.floatValue = 0.2f;
	maxSeed.floatValue = 1.0f;

	minSeed3.floatValue = -2.22554f;
	maxSeed3.floatValue = -1.0f;
	
	minSeed2.uint8Value = 0;
	maxSeed2.uint8Value = 1;

	Gene *gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed, maxSeed);
	gene->enableBounding(true);
	gene->setBounds(minSeed, maxSeed);
	genotype->push_back(gene);

	gene = new Gene(FLOAT);
	gene->setSeedRange(minSeed3, maxSeed3);
	gene->enableBounding(true);
	gene->setBounds(minSeed3, maxSeed3);
	genotype->push_back(gene);

	gene = new Gene(UINT8);
	gene->setSeedRange(minSeed2, maxSeed2);
	gene->enableBounding(true);
	gene->setBounds(minSeed2, maxSeed2);
	genotype->push_back(gene);

	return genotype;
}

double fitnessFunction3(Chromosome * chromosome) {
	//MODE_MAXIMIZATION	

	float x1 = (*chromosome->getGenes())[0]->getValue().floatValue;	
	float x2 = (*chromosome->getGenes())[1]->getValue().floatValue;
	float y = (float)(*chromosome->getGenes())[2]->getValue().uint8Value;

	//We have to check if the float/double values are valid
	if (isnan(x1) || isinf(x1))
		return -INFINITY;
	if (isnan(x2) || isinf(x2))
		return -INFINITY;

	double fitness = -(-0.7f*y + 5.0f*pow(x1 - 0.5f, 2.0f) + 0.8f);
	double punish = 0;

	bool violated = false;
	float violations[3] = { 0, 0, 0};	

	if (-exp(x1 - 0.2f) - x2 > 0) {
		violated = true;
		violations[0] = -exp(x1 - 0.2f) - x2;
	}
	if (x2 + 1.1f*y > -1.0f) {
		violated = true;
		violations[1] = -1.0f - (x2 + 1.1f*y);
	}
	if (x1 - 1.2f*y > 0.2f) {
		violated = true;
		violations[2] = 0.2f - (x1 - 1.2f*y);
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <3; i++)
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}
	
	return fitness + punish;
}

vector<Gene*> *getGenotype2() {
	//MODE_MAXIMIZATION
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

double fitnessFunction2(Chromosome * chromosome) {
	//MODE_MAXIMIZATION	

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
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}

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

double fitnessFunction2d(Chromosome * chromosome) {
	//MODE_MAXIMIZATION

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
			punish += abs(violations[i]);
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
	//MODE_MAXIMIZATION

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
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}

	return fitness + punish;
}