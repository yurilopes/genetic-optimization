#pragma once

#include "mutationvectorized.h"
#include "chromosome.h"
#include <random>
#include <math.h>

class MutationGaussian : public MutationVectorized {
	protected:
		double									mean;
		double									stddev;
		std::normal_distribution<double>		*distribution = NULL;
		std::mt19937							*randomGen = NULL;
		void									initialize(double m, double sdev);

		void									refreshDistribution();

	public:
		void									mutate(Chromosome *chromosome);
		void									mutate(Gene *gene);
		double									getMean();
		void									setMean(double m);
		double									getStdDev();
		void									setStdDev(double s);

		MutationGaussian(double m, double sdev);		
		~MutationGaussian();
		
};

inline void MutationGaussian::initialize(double m, double sdev)
{
	mean = m;
	stddev = sdev;
	randomGen = new std::mt19937(static_cast<unsigned int>(std::time(0)));
	distribution = new std::normal_distribution<double>(mean, stddev);
}

inline void MutationGaussian::refreshDistribution()
{
	if (distribution == NULL)
		return;

	delete distribution;

	distribution = new std::normal_distribution<double>(mean, stddev);
}

inline void MutationGaussian::mutate(Chromosome* chromosome)
{

	std::vector<Gene *> *genes = chromosome->getGenes();
	size_t i = 0;
	
	for (std::vector<Gene *>::iterator it = genes->begin(); it != genes->end(); it++) {		

		mutate(*it);

		i++;
	}
}

inline void MutationGaussian::mutate(Gene* gene) {
	double randomNum = (*distribution)(*randomGen);
	GeneValue newValue;

	switch (gene->getDataType())
	{
		case (FLOAT): {
			newValue.floatValue = gene->getValue().floatValue + (float)randomNum;
			gene->setValue(newValue);
			break;
		}
		case (DOUBLE): {
			newValue.doubleValue = gene->getValue().doubleValue + randomNum;
			gene->setValue(newValue);
			break;
		}
		case (INT8): {
			newValue.int8Value = gene->getValue().int8Value + (int8_t)round(randomNum);
			gene->setValue(newValue);
			break;
		}
		case (UINT8): {
			newValue.uint8Value = (uint8_t)(gene->getValue().int8Value + (int8_t)round(randomNum));
			gene->setValue(newValue);
			break;
		}
		case (INT16): {
			newValue.int16Value = gene->getValue().int16Value + (int16_t)round(randomNum);
			gene->setValue(newValue);
			break;
		}
		case (UINT16): {
			newValue.uint16Value = (uint16_t)(gene->getValue().int16Value + (int16_t)round(randomNum));
			gene->setValue(newValue);
			break;
		}
		case (INT32): {
			newValue.int32Value = gene->getValue().int32Value + (int32_t)round(randomNum);
			gene->setValue(newValue);
			break;
		}
		case (UINT32): {
			newValue.uint32Value = (uint32_t)(gene->getValue().int32Value + (int32_t)round(randomNum));
			gene->setValue(newValue);
			break;
		}
		case (INT64): {
			newValue.int64Value = gene->getValue().int64Value + (int64_t)round(randomNum);
			gene->setValue(newValue);
			break;
		}
		case (UINT64): {
			newValue.uint64Value = (uint64_t)(gene->getValue().int64Value + (int64_t)round(randomNum));
			gene->setValue(newValue);
			break;
		}
	}
}

inline double MutationGaussian::getMean()
{
	return mean;
}

inline void MutationGaussian::setMean(double m)
{
	mean = m;
	refreshDistribution();
}

inline double MutationGaussian::getStdDev()
{
	return stddev;
}

inline void MutationGaussian::setStdDev(double s)
{
	stddev = s;
	refreshDistribution();
}

inline MutationGaussian::MutationGaussian(double m, double sdev)
{
	initialize(m, sdev);
}

inline MutationGaussian::~MutationGaussian()
{
	if (distribution)
		delete distribution;
	if (randomGen)
		delete randomGen;
}
