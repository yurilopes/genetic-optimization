#pragma once

#include <random>

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
	DOUBLE
};


/*
This union holds the value within a Gene

It doesn't matter what data type you want to use ou read the gene value as.
This means you can create and work with a float gene but at any point read its value as a int16_t as well
*/
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
	GeneValue				lowerBound;
	GeneValue				upperBound;	
	bool					boundingEnabled = false;

public:
	Gene(GeneDataType dataTy);
	Gene(Gene * original);
	GeneValue				getValue();
	GeneValue				*getValuePointer(); //Mostly used for bitwise genetic operators
	void					setValue(GeneValue val);
	GeneDataType			getDataType();

	void					setSeedRange(GeneValue min, GeneValue max);
	void					setBounds(GeneValue lower, GeneValue upper);
	void					refreshValue(); //Useful when resetting the boundings of a gene
	void					enableBounding(bool enable);

	bool					isBoundingEnabled();
	GeneValue				getMinimumSeed();
	GeneValue				getMaximumSeed();
	GeneValue				getLowerBound();
	GeneValue				getUpperBound();	

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

inline GeneValue * Gene::getValuePointer()
{
	return &value;
}

inline void Gene::setValue(GeneValue val)
{
	value = val;
	if (!boundingEnabled)
		return;

	//If bounding is enabled we should clip the value if it goes out of bounds

	switch (dataType)
	{
		case (FLOAT): {
			if (value.floatValue < lowerBound.floatValue)
				value = lowerBound;
			if (value.floatValue > upperBound.floatValue)
				value = upperBound;
			break;
		}
		case (DOUBLE): {
			if (value.doubleValue < lowerBound.doubleValue)
				value = lowerBound;
			if (value.doubleValue > upperBound.doubleValue)
				value = upperBound;
			break;
		}
		case (INT8): {
			if (value.int8Value < lowerBound.int8Value)
				value = lowerBound;
			if (value.int8Value > upperBound.int8Value)
				value = upperBound;
			break;
		}
		case (UINT8): {
			if (value.uint8Value < lowerBound.uint8Value)
				value = lowerBound;
			if (value.uint8Value > upperBound.uint8Value)
				value = upperBound;
			break;
		}
		case (INT16): {
			if (value.int16Value < lowerBound.int16Value)
				value = lowerBound;
			if (value.int16Value > upperBound.int16Value)
				value = upperBound;
			break;
		}
		case (UINT16): {
			if (value.uint16Value < lowerBound.uint16Value)
				value = lowerBound;
			if (value.uint16Value > upperBound.uint16Value)
				value = upperBound;
			break;
		}
		case (INT32): {
			if (value.int32Value < lowerBound.int32Value)
				value = lowerBound;
			if (value.int32Value > upperBound.int32Value)
				value = upperBound;
			break;
		}
		case (UINT32): {
			if (value.uint32Value < lowerBound.uint32Value)
				value = lowerBound;
			if (value.uint32Value > upperBound.uint32Value)
				value = upperBound;
			break;
		}
		case (INT64): {
			if (value.int64Value < lowerBound.int64Value)
				value = lowerBound;
			if (value.int64Value > upperBound.int64Value)
				value = upperBound;
			break;
		}
		case (UINT64): {
			if (value.uint64Value < lowerBound.uint64Value)
				value = lowerBound;
			if (value.uint64Value > upperBound.uint64Value)
				value = upperBound;
			break;
		}
	}

}

inline GeneDataType Gene::getDataType()
{
	return dataType;
}

inline void Gene::setSeedRange(GeneValue min, GeneValue max)
{
	minSeed = min;
	maxSeed = max;
}

inline void Gene::setBounds(GeneValue lower, GeneValue upper)
{
	lowerBound = lower;
	upperBound = upper;
}

inline void Gene::refreshValue()
{
	setValue(value);
}

inline void Gene::enableBounding(bool enable)
{
	boundingEnabled = enable;
}

inline bool Gene::isBoundingEnabled()
{
	return boundingEnabled;
}

inline GeneValue Gene::getMinimumSeed()
{
	return minSeed;
}

inline GeneValue Gene::getMaximumSeed()
{
	return maxSeed;
}

inline GeneValue Gene::getLowerBound()
{
	return lowerBound;
}

inline GeneValue Gene::getUpperBound()
{
	return upperBound;
}

inline void Gene::initialize()
{
	switch (dataType)
	{
		case (FLOAT): {
			uniform_real_distribution<float> randomF(minSeed.floatValue, maxSeed.floatValue);
			value.floatValue = randomF(genGA);
			break;
		}
		case (DOUBLE): {
			uniform_real_distribution<double> randomD(minSeed.doubleValue, maxSeed.doubleValue);
			value.doubleValue = randomD(genGA);
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
	}
}
