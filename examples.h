#pragma once

#define	PUNISHMENT_FACTOR		-1e3f

#include "genetic-algorithm.h"
#include <algorithm>

#define EX8_M	6
#define EX8_N	5
#define EX8_H	6000.0

double S8[EX8_N][EX8_M] = {
	{ 7.9, 2.0, 5.2, 4.9, 6.1, 4.2 },
	{ 0.7, 0.8, 0.9, 3.4, 2.1, 2.5 },
	{ 0.7, 2.6, 1.6, 3.6, 3.2, 2.9 },
	{ 4.7, 2.3, 1.6, 2.7, 1.2, 2.5 },
	{ 1.2, 3.6, 2.4, 4.5, 1.6, 2.1 }
};

double t8[EX8_N][EX8_M] = {
	{ 6.4, 4.7, 8.3, 3.9, 2.1, 1.2 },
	{ 6.8, 6.4, 6.5, 4.4, 2.3, 3.2 },
	{ 1.0, 6.3, 5.4, 11.9, 5.7, 6.2 },
	{ 3.2, 3.0, 3.5, 3.3, 2.8, 3.4 },
	{ 2.1, 2.5, 4.2, 3.6, 3.7, 2.2 }
};

double Q8[EX8_N] = { 250000.0, 150000.0, 180000.0, 160000.0, 120000.0 };

double TLiMIN8(int i, int M, double Nju) {
	bool first = true;

	double maxValue;

	for (int j = 0; j < M; j++) {
		if (first) {
			first = false;
			maxValue = t8[i][j]/Nju;
		}
		if(t8[i][j] / Nju> maxValue)
			maxValue = t8[i][j] / Nju;
	}

	return maxValue;
}

double TLiMAX8(int i, int M) {
	bool first = true;

	double maxValue;

	for (int j = 0; j < M; j++) {
		if (first) {
			first = false;
			maxValue = t8[i][j];
		}
		if (t8[i][j]> maxValue)
			maxValue = t8[i][j];
	}

	return maxValue;
}

double BiMIN8(int i, int M) {
	return (Q8[i] / EX8_H) * TLiMIN8(i, M, 4.0);
}

double BiMAX8(int i, int M) {
	double minVS;
	bool first = true;

	for (int j = 0; j < M; j++) {
		if (first) {
			first = false;
			minVS = 3000.0 / S8[i][j];
		}
		if (3000.0 / S8[i][j] < minVS)
			minVS = 3000.0 / S8[i][j];
	}

	return std::min(Q8[i], minVS);

}

vector<Gene*> *getGenotype8() {
	GeneValue minSeed, maxSeed;
	vector<Gene *> *genotype = new vector<Gene *>();

	Gene * gene;

	minSeed.uint8Value = 1;
	maxSeed.uint8Value = 4;

	//N1..N6
	for (int i = 0; i < EX8_M; i++) {
		/*
		if (i == 0) {
			minSeed.uint8Value = 2;
			maxSeed.uint8Value = 2;
		}
		if (i == 1) {
			minSeed.uint8Value = 2;
			maxSeed.uint8Value = 2;
		}
		if (i == 2) {
			minSeed.uint8Value = 3;
			maxSeed.uint8Value = 3;
		}
		if (i == 3) {
			minSeed.uint8Value = 2;
			maxSeed.uint8Value = 2;
		}
		if (i == 4) {
			minSeed.uint8Value = 1;
			maxSeed.uint8Value = 1;
		}
		if (i == 5) {
			minSeed.uint8Value = 1;
			maxSeed.uint8Value = 1;
		}
		*/

		gene = new Gene(UINT8);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	minSeed.floatValue = 300.0f;
	maxSeed.floatValue = 3000.0f;

	//V1..V6
	for (int i = 0; i < EX8_M; i++) {
		minSeed.floatValue = 300.0f;
		maxSeed.floatValue = 3000.0f;
		if (i == 0) {
			minSeed.floatValue = 3000.0;
			maxSeed.floatValue = 3000.0;
		}
		if (i == 1) {
			minSeed.floatValue = 1892.0;
			maxSeed.floatValue = 1892.0;
		}
		if (i == 2) {
			minSeed.floatValue = 1975.0;
			maxSeed.floatValue = 1975.0;
		}
		
		if (i == 3) {
			minSeed.floatValue = 2619.0;
			maxSeed.floatValue = 2619.0;
		}
		if (i == 4) {
			minSeed.floatValue = 2328.0;
			maxSeed.floatValue = 2328.0;
		}
		if (i == 5) {
			minSeed.floatValue = 2110.0;
			maxSeed.floatValue = 2110.0;
		}
		

		gene = new Gene(FLOAT);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	//B1..B5
	for (int i = 0; i < EX8_N; i++) {
		minSeed.floatValue = (float)BiMIN8(i, EX8_M);
		maxSeed.floatValue = (float)BiMAX8(i, EX8_M);

		gene = new Gene(FLOAT);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	//TL1..TL5
	for (int i = 0; i < EX8_N; i++) {
		minSeed.floatValue = (float)TLiMIN8(i, EX8_M, 4.0);
		maxSeed.floatValue = (float)TLiMAX8(i, EX8_M);

		gene = new Gene(FLOAT);
		gene->setSeedRange(minSeed, maxSeed);
		gene->enableBounding(true);
		gene->setBounds(minSeed, maxSeed);
		genotype->push_back(gene);
	}

	return genotype;
}

double fitnessFunction8(Chromosome * chromosome) {

	double N[EX8_M];
	double V[EX8_M];
	double B[EX8_N];
	double TL[EX8_N];

	N[0] = (double)(*chromosome->getGenes())[0]->getValue().uint8Value;
	N[1] = (double)(*chromosome->getGenes())[1]->getValue().uint8Value;
	N[2] = (double)(*chromosome->getGenes())[2]->getValue().uint8Value;
	N[3] = (double)(*chromosome->getGenes())[3]->getValue().uint8Value;
	N[4] = (double)(*chromosome->getGenes())[4]->getValue().uint8Value;
	N[5] = (double)(*chromosome->getGenes())[5]->getValue().uint8Value;

	V[0] = (double)(*chromosome->getGenes())[6]->getValue().floatValue;
	V[1] = (double)(*chromosome->getGenes())[7]->getValue().floatValue;
	V[2] = (double)(*chromosome->getGenes())[8]->getValue().floatValue;
	V[3] = (double)(*chromosome->getGenes())[9]->getValue().floatValue;
	V[4] = (double)(*chromosome->getGenes())[10]->getValue().floatValue;
	V[5] = (double)(*chromosome->getGenes())[11]->getValue().floatValue;

	B[0] = (double)(*chromosome->getGenes())[12]->getValue().floatValue;
	B[1] = (double)(*chromosome->getGenes())[13]->getValue().floatValue;
	B[2] = (double)(*chromosome->getGenes())[14]->getValue().floatValue;
	B[3] = (double)(*chromosome->getGenes())[15]->getValue().floatValue;
	B[4] = (double)(*chromosome->getGenes())[16]->getValue().floatValue;

	TL[0] = (double)(*chromosome->getGenes())[17]->getValue().floatValue;
	TL[1] = (double)(*chromosome->getGenes())[18]->getValue().floatValue;
	TL[2] = (double)(*chromosome->getGenes())[19]->getValue().floatValue;
	TL[3] = (double)(*chromosome->getGenes())[20]->getValue().floatValue;
	TL[4] = (double)(*chromosome->getGenes())[21]->getValue().floatValue;

	double fitness = 0;
	for (int j = 0; j < EX8_M; j++)
		fitness -= 250.0 * N[j] * pow(V[j], 0.6);

	double punish = 0;

	bool violated = false;
	double violations[61];

	for (int i = 0; i < 61; i++)
		violations[i] = 0;

	int iv = 0;
	double sum = 0;
	for (int i = 0; i < EX8_N; i++) {
		sum += (Q8[i] * TL[i]) / B[i];
	}

	if (sum > EX8_H) {
		violated = true;
		violations[iv] = EX8_H - sum;
	}

	iv++;

	for (int i = 0; i < EX8_N; i++) {
		for (int j = 0; j < EX8_M; j++) {
			if (V[j] < S8[i][j] * B[i]) {
				violated = true;
				violations[iv] = (S8[i][j] * B[i]) - (V[j]);
			}
			iv++;
		}
	}

	for (int i = 0; i < EX8_N; i++) {
		for (int j = 0; j < EX8_M; j++) {
			if (N[j] * TL[i] < t8[i][j]) {
				violated = true;
				violations[iv] = t8[i][j] - (N[j] * TL[i]);
			}
			iv++;
		}
	}

	if (violated) {
		punish = 0;
		for (int i = 0; i <61; i++)
			punish += abs(violations[i]);
		punish *= PUNISHMENT_FACTOR;
	}


	return fitness + punish;
}

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
	double sum = 0;
	for (int i = 0; i < 2; i++) {
		sum += (Q[i] * TL[i]) / B[i];
	}

	if (sum > 6000.0) {
		violated = true;
		violations[iv] = 6000.0 - sum;
	}

	iv++;

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
	if (std::isnan(x1) || std::isinf(x1))
		return -INFINITY;
	if (std::isnan(x2) || std::isinf(x2))
		return -INFINITY;
	if (std::isnan(x3) || std::isinf(x3))
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
	if (std::isnan(x1) || std::isinf(x1))
		return -INFINITY;
	if (std::isnan(x2) || std::isinf(x2))
		return -INFINITY;
	if (std::isnan(x3) || std::isinf(x3))
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
	if (std::isnan(v1) || std::isinf(v1))
		return INFINITY;
	if (std::isnan(v2) || std::isinf(v2))
		return INFINITY;
	if (std::isnan(x1) || std::isinf(x1))
		return INFINITY;
	if (std::isnan(x2) || std::isinf(x2))
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
	if (std::isnan(x1) || std::isinf(x1))
		return -INFINITY;
	if (std::isnan(x2) || std::isinf(x2))
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
	if (std::isnan(x1) || std::isinf(x1))
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
	if (std::isnan(x1) || std::isinf(x1))
		return -INFINITY;
	if (std::isnan(x2) || std::isinf(x2))
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
	if (std::isnan(x) || std::isinf(x))
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