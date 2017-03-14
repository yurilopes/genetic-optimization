#pragma once

#include <vector>
#include <climits>
#include <iostream>

extern mt19937 genGA;

enum GeneDataType {
	UINT8,
	INT8,
	UINT16,
	INT16,
	UINT32,
	INT32,
	UINT64,
	INT64,
	FLOAT,
	DOUBLE,
	CUSTOM
};

union GeneValue {
	float		floatValue;
	double		doubleValue;
	int8_t		int8Value;
	int16_t		int16Value;
	int32_t		int32Value;
	int64_t		int64Value;
	uint8_t		uint8Value;
	uint16_t	uint16Value;
	uint32_t	uint32Value;
	uint64_t	uint64Value;
};

class Gene {

	protected:
		GeneValue				value;
		GeneDataType			dataType;		
		GeneValue				minSeed;
		GeneValue				maxSeed;
		
		//Properties for CUSTOM data type
		uint64_t				lowerBound, upperBound;

	public:
		Gene(GeneDataType dataTy);
		Gene(Gene * original);		
		GeneValue				getValue();
		void					setValue(GeneValue val);
		GeneDataType			getDataType();		

		void					setSeedRange(uint64_t lower, uint64_t upper);
		void					setSeedRange(int64_t lower, int64_t upper);
		void					setSeedRange(uint32_t lower, uint32_t upper);
		void					setSeedRange(int32_t lower, int32_t upper);
		void					setSeedRange(uint16_t lower, uint16_t upper);
		void					setSeedRange(int16_t lower, int16_t upper);
		void					setSeedRange(uint8_t lower, uint8_t upper);
		void					setSeedRange(int8_t lower, int8_t upper);
		void					setSeedRange(float lower, float upper);
		void					setSeedRange(double lower, double upper);	

		GeneValue				getMinimumSeed();
		GeneValue				getMaximumSeed();	

		void					initialize();
};

Gene::Gene(GeneDataType dataTy) {
	dataType = dataTy;

	//Either double or uint64_t will be the biggest data field
	if (sizeof(double) > sizeof(uint64_t))
		value.doubleValue = 0;
	else
		value.uint64Value = 0;	
}

inline Gene::Gene(Gene *original)
{
	dataType = original->getDataType();	

	//Either double or uint64_t will be the biggest data field
	minSeed = original->getMinimumSeed();
	maxSeed = original->getMaximumSeed();
	value = original->getValue();
}

inline GeneValue Gene::getValue()
{
	return value;
}

inline void Gene::setValue(GeneValue val)
{
	value = val;
}

inline GeneDataType Gene::getDataType()
{
	return dataType;
}

inline void Gene::setSeedRange(uint64_t lower, uint64_t upper)
{
	minSeed.uint64Value = lower;
	maxSeed.uint64Value = upper;
}

inline void Gene::setSeedRange(int64_t lower, int64_t upper)
{
	minSeed.int64Value = lower;
	maxSeed.int64Value = upper;
}

inline void Gene::setSeedRange(uint32_t lower, uint32_t upper)
{
	minSeed.uint32Value = lower;
	maxSeed.uint32Value = upper;
}

inline void Gene::setSeedRange(int32_t lower, int32_t upper)
{
	minSeed.int32Value = lower;
	maxSeed.int32Value = upper;
}

inline void Gene::setSeedRange(uint16_t lower, uint16_t upper)
{
	minSeed.uint16Value = lower;
	maxSeed.uint16Value = upper;
}

inline void Gene::setSeedRange(int16_t lower, int16_t upper)
{
	minSeed.int16Value = lower;
	maxSeed.int16Value = upper;
}

inline void Gene::setSeedRange(uint8_t lower, uint8_t upper)
{
	minSeed.uint8Value = lower;
	maxSeed.uint8Value = upper;
}

inline void Gene::setSeedRange(int8_t lower, int8_t upper)
{
	minSeed.int8Value = lower;
	maxSeed.int8Value = upper;
}

inline void Gene::setSeedRange(float lower, float upper)
{
	minSeed.floatValue = lower;
	maxSeed.floatValue = upper;
}

inline void Gene::setSeedRange(double lower, double upper)
{
	minSeed.doubleValue = lower;
	minSeed.doubleValue = upper;
}

inline GeneValue Gene::getMinimumSeed()
{
	return minSeed;
}

inline GeneValue Gene::getMaximumSeed()
{
	return maxSeed;
}

inline void Gene::initialize()
{
	switch (dataType)
	{
		case (FLOAT): { 
			uniform_real_distribution<float> randomF(minSeed.floatValue, maxSeed.floatValue);
			GeneValue val;
			val.floatValue = randomF(genGA);
			value = val;
			break;
		}
		case (DOUBLE): {
			uniform_real_distribution<double> randomD(minSeed.doubleValue, maxSeed.doubleValue);
			GeneValue val;
			val.doubleValue = randomD(genGA);
			value = val;
			break;
		}
		case (INT8): {
			uniform_int_distribution<int16_t> randomI16(minSeed.int8Value, maxSeed.int8Value);
			#pragma warning(push)
			#pragma warning(disable : 4244)
			/*
			Should never cause loss of data, since seed values are within the 8 bit range
			*/
			GeneValue val;
			val.int8Value = randomI16(genGA);
			value = val;
			#pragma warning(pop)
			break;
		}
		case (UINT8): {
			/*
			Can't have a uniform_int_distribution<uint8_t> because:

			It seems that this library implementaion behaviour does not violate the ISO C++ standard, because (26.5.1.1)
			"effect of instantiating a template ... that has a template type parameter named IntType __is undefined unless__ the
			corresponding template argument is cv-unqualified
			and __is one of short, int, long, long long, unsigned short, unsigned int, unsigned long, or unsigned long long__",
			char is not listed.

			More info at: http://stackoverflow.com/questions/31460733/why-arent-stduniform-int-distributionuint8-t-and-stduniform-int-distri
			*/
			uniform_int_distribution<uint16_t> randomI16(minSeed.uint8Value, maxSeed.uint8Value);
			#pragma warning(push)
			#pragma warning(disable : 4244)
			/*
			Should never cause loss of data, since seed values are within the 8 bit range
			*/
			GeneValue val;
			val.uint8Value = randomI16(genGA);
			value = val;
			#pragma warning(pop)
			break;
		}
		case (INT16): {
			uniform_int_distribution<int16_t> randomI16(minSeed.int16Value, maxSeed.int16Value);
			GeneValue val;
			val.int16Value = randomI16(genGA);
			value = val;
			break;
		}
		case (UINT16): {
			uniform_int_distribution<uint16_t> randomI16(minSeed.uint16Value, maxSeed.uint16Value);
			GeneValue val;
			val.uint16Value = randomI16(genGA);
			value = val;
			break;
		}
		case (INT32): {
			uniform_int_distribution<int32_t> randomI32(minSeed.int32Value, maxSeed.int32Value);
			GeneValue val;
			val.int32Value = randomI32(genGA);
			value = val;
			break;
		}
		case (UINT32): {
			uniform_int_distribution<uint32_t> randomI32(minSeed.uint32Value, maxSeed.uint32Value);
			GeneValue val;
			val.uint32Value = randomI32(genGA);
			value = val;
			break;
		}
		case (INT64): {
			uniform_int_distribution<int64_t> randomI64(minSeed.int64Value, maxSeed.int64Value);
			GeneValue val;
			val.int64Value = randomI64(genGA);
			value = val;
			break;
		}
		case (UINT64): {
			uniform_int_distribution<uint64_t> randomI64(minSeed.uint64Value, maxSeed.uint64Value);
			GeneValue val;
			val.uint64Value = randomI64(genGA);
			value = val;
			break;
		}
		default: {
			//TODO: Set value for CUSTOM data type
			break;
		}
	}
}
