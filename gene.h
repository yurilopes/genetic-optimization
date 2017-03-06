#pragma once

#include <vector>
#include <climits>

#if CHAR_BIT != 8
	#pragma message("CHAR_BIT is not 8 bits long. This could lead to erratic behaviour and loss of information.")
#endif // CHAR_BIT != 8


enum GeneDataType {
	INT8,
	INT16,
	INT32,
	INT64,
	FLOAT,
	DOUBLE,
	CUSTOM
};

class Gene {

	protected:
		union SeedValue {
			float		floatValue;
			double		doubleValue;
			uint64_t	uintValue;
		};

		vector<unsigned char>	*geneBits = NULL;
		GeneDataType			dataType;		
		SeedValue				minSeed, maxSeed;
		
		//Properties for CUSTOM data type
		uint64_t				lowerBound, upperBound;

	public:
		Gene(GeneDataType dataTy);
		~Gene();
		GeneDataType	getDataType();
		uint8_t			getValueUInt8();
		int8_t			getValueInt8();
		uint16_t		getValueUInt16();
		int16_t			getValueInt16();
		uint32_t		getValueUInt32();
		int32_t			getValueInt32();
		uint64_t		getValueUInt64();
		int64_t			getValueInt64();
		float			getValueFloat();
		double			getValueDouble();
		void			setSeedRange(uint64_t lower, uint64_t upper);
		void			setSeedRange(int64_t lower, int64_t upper);
		void			setSeedRange(float lower, float upper);
		void			setSeedRange(double lower, double upper);

};

Gene::Gene(GeneDataType dataTy) {
	dataType = dataTy;

	/*
	CUSTOM data type is not previously allocated since it requires upper and lower bounds to be defined
	*/
	
	switch (dataType)
	{
		case INT8: 
			geneBits = new vector<unsigned char>(sizeof(uint8_t) * CHAR_BIT, 0);
			break;
		case INT16:
			geneBits = new vector<unsigned char>(sizeof(uint16_t) * CHAR_BIT, 0);
			break;
		case INT32:
			geneBits = new vector<unsigned char>(sizeof(uint32_t) * CHAR_BIT, 0);
			break;
		case INT64:
			geneBits = new vector<unsigned char>(sizeof(uint64_t) * CHAR_BIT, 0);
			break;
		case FLOAT:
			geneBits = new vector<unsigned char>(sizeof(float) * CHAR_BIT, 0);
			if (sizeof(float) != sizeof(uint32_t))
				cerr << "float is not " << sizeof(uint32_t) << " bytes long  (same as uint32_t). This could lead to erratic behaviour and loss of information." << endl;
			break;
		case DOUBLE:
			geneBits = new vector<unsigned char>(sizeof(double) * CHAR_BIT, 0);
			if (sizeof(double) != sizeof(uint64_t))
				cerr << "double is not " << sizeof(uint64_t) << " bytes long  (same as uint64_t). This could lead to erratic behaviour and loss of information." << endl;
			break;
		default:
			break;
	}
}

Gene::~Gene() {
	if (geneBits)
		delete geneBits;
}

inline GeneDataType Gene::getDataType()
{
	return dataType;
}

inline uint8_t Gene::getValueUInt8()
{
	uint8_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (i == sizeof(uint8_t)*CHAR_BIT)
			break;
		if (bv)
			result = result | (1 << i);
		i++;
	}
	return result;
}

inline int8_t Gene::getValueInt8()
{
	uint8_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (i == sizeof(uint8_t)*CHAR_BIT)
			break;
		if (bv)
			result = result | (1 << i);
		i++;
	}
	return (int8_t)result;
}

inline uint16_t Gene::getValueUInt16()
{
	uint16_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (i == sizeof(uint16_t)*CHAR_BIT)
			break;
		if (bv)
			result = result | (1 << i);
		i++;
	}
	return result;
}

inline int16_t Gene::getValueInt16()
{
	uint16_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (i == sizeof(uint16_t)*CHAR_BIT)
			break;
		if (bv)
			result = result | (1 << i);
		i++;
	}
	return (int16_t)result;
}

inline uint32_t Gene::getValueUInt32()
{
	uint32_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (i == sizeof(uint32_t)*CHAR_BIT)
			break;
		if (bv)
			result = result | (1 << i);
		i++;
	}
	return result;
}

inline int32_t Gene::getValueInt32()
{
	uint32_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (i == sizeof(uint32_t)*CHAR_BIT)
			break;
		if (bv)
			result = result | (1 << i);
		i++;
	}
	return (int32_t)result;
}

inline uint64_t Gene::getValueUInt64()
{
	uint64_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (bv)
			result = result | (1ui64 << i);
		i++;
	}
	return result;
}

inline int64_t Gene::getValueInt64()
{
	uint64_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (bv)
			result = result | (1ui64 << i);
		i++;
	}
	return (int64_t)result;
}

inline float Gene::getValueFloat()
{
	if (sizeof(float) != sizeof(uint32_t))
		cerr << "float is not " << sizeof(uint32_t) << " bytes long  (same as uint32_t). This could lead to erratic behaviour and loss of information." << endl;

	uint32_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (i == sizeof(uint32_t)*CHAR_BIT)
			break;
		if (bv)
			result = result | (1 << i);
		i++;
	}
	
	/*
	Memory cast from pointer to uint32_t to float pointer at same memory position
	This makes the code interpret the bytes at &result as a float
	*/
	return *(float *)(&result); 
}

inline double Gene::getValueDouble()
{
	if (sizeof(double) != sizeof(uint64_t))
		cerr << "double is not " << sizeof(uint64_t) << " bytes long (same as uint64_t). This could lead to erratic behaviour and loss of information." << endl;

	uint64_t result = 0;
	unsigned int i = 0;
	for (unsigned char bv : *geneBits) {
		if (bv)
			result = result | (1ui64 << i);
		i++;
	}

	/*
	Memory cast from pointer to uint64_t to double pointer at same memory position
	This makes the code interpret the bytes at &result as a double
	*/
	return *(double *)(&result);
}

inline void Gene::setSeedRange(uint64_t lower, uint64_t upper)
{
	minSeed.uintValue = lower;
	maxSeed.uintValue = upper;
}

inline void Gene::setSeedRange(int64_t lower, int64_t upper)
{
	minSeed.uintValue = (uint64_t)lower;
	maxSeed.uintValue = (uint64_t)upper;
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
