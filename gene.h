#pragma once

#include <vector>
#include <climits>
#include <iostream>

#define ERR_DOUBLE_SIZEMISMATCH "double is not " << sizeof(uint64_t) << " bytes long  (same as uint64_t). This could lead to erratic behaviour and loss of information."
#define ERR_FLOAT_SIZEMISMATCH "float is not " << sizeof(uint32_t) << " bytes long  (same as uint32_t). This could lead to erratic behaviour and loss of information."

#if CHAR_BIT != 8
	#pragma message("CHAR_BIT is not 8 bits long. This could lead to erratic behaviour and loss of information.")
#endif // CHAR_BIT != 8


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
		vector<unsigned char>	*geneBits = NULL;
		GeneDataType			dataType;		
		GeneValue				minSeed;
		GeneValue				maxSeed;
		
		//Properties for CUSTOM data type
		uint64_t				lowerBound, upperBound;

	public:
		Gene(GeneDataType dataTy);
		Gene(Gene * original);
		~Gene();
		vector<unsigned char> *getGeneBits();
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

		/*
		There may be a better way to do these below
		However, the day I wrote these I was sick and not thinking straight
		Even if these functions take some, maybe, unnecessary space, they work as intended
		*/
		void			setValueInt8(int8_t value); //TODO
		void			setValueUInt8(uint8_t value);
		void			setValueInt16(int16_t value);
		void			setValueUInt16(uint16_t value);
		void			setValueInt32(int32_t value);
		void			setValueUInt32(uint32_t value);
		void			setValueInt64(int64_t value);
		void			setValueUInt64(uint64_t value);
		void			setValueFloat(float value);
		void			setValueDouble(double value);
		void			setSeedRange(uint64_t lower, uint64_t upper);
		void			setSeedRange(int64_t lower, int64_t upper);
		void			setSeedRange(uint32_t lower, uint32_t upper);
		void			setSeedRange(int32_t lower, int32_t upper);
		void			setSeedRange(uint16_t lower, uint16_t upper);
		void			setSeedRange(int16_t lower, int16_t upper);
		void			setSeedRange(uint8_t lower, uint8_t upper);
		void			setSeedRange(int8_t lower, int8_t upper);
		void			setSeedRange(float lower, float upper);
		void			setSeedRange(double lower, double upper);	

		GeneValue		getMinimumSeed();
		GeneValue		getMaximumSeed();

		void			printBits();
};

Gene::Gene(GeneDataType dataTy) {
	dataType = dataTy;

	/*
	CUSTOM data type is not previously allocated since it requires upper and lower bounds to be defined
	*/
	
	switch (dataType)
	{
		case INT8: 
		case UINT8:
			geneBits = new vector<unsigned char>(sizeof(uint8_t) * CHAR_BIT, 0);
			break;
		case INT16:
		case UINT16:
			geneBits = new vector<unsigned char>(sizeof(uint16_t) * CHAR_BIT, 0);
			break;
		case INT32:
		case UINT32:
			geneBits = new vector<unsigned char>(sizeof(uint32_t) * CHAR_BIT, 0);
			break;
		case INT64:
		case UINT64:
			geneBits = new vector<unsigned char>(sizeof(uint64_t) * CHAR_BIT, 0);
			break;
		case FLOAT:			
			if (sizeof(float) != sizeof(uint32_t))
				cerr << ERR_FLOAT_SIZEMISMATCH << endl;
			geneBits = new vector<unsigned char>(sizeof(float) * CHAR_BIT, 0);
			break;
		case DOUBLE:			
			if (sizeof(double) != sizeof(uint64_t))
				cerr << ERR_DOUBLE_SIZEMISMATCH << endl;
			geneBits = new vector<unsigned char>(sizeof(double) * CHAR_BIT, 0);
			break;
		default: //TODO
			break;
	}
}

inline Gene::Gene(Gene *original) //Construction from gene model
{
	dataType = original->getDataType();	

	switch (dataType)
	{
		case INT8:
		case UINT8:
			minSeed.uint8Value = original->getMinimumSeed().uint8Value;
			maxSeed.uint8Value = original->getMaximumSeed().uint8Value;
			geneBits = new vector<unsigned char>(*original->getGeneBits());
			break;
		case INT16:
		case UINT16:
			minSeed.uint16Value = original->getMinimumSeed().uint16Value;
			maxSeed.uint16Value = original->getMaximumSeed().uint16Value;
			geneBits = new vector<unsigned char>(*original->getGeneBits());
			break;
		case INT32:
		case UINT32:
			minSeed.uint32Value = original->getMinimumSeed().uint32Value;
			maxSeed.uint32Value = original->getMaximumSeed().uint32Value;
			geneBits = new vector<unsigned char>(*original->getGeneBits());
			break;
		case INT64:
		case UINT64:
			minSeed.uint64Value = original->getMinimumSeed().uint64Value;
			maxSeed.uint64Value = original->getMaximumSeed().uint64Value;
			geneBits = new vector<unsigned char>(*original->getGeneBits());
			break;
		case FLOAT:
			minSeed.floatValue = original->getMinimumSeed().floatValue;
			maxSeed.floatValue = original->getMaximumSeed().floatValue;
			geneBits = new vector<unsigned char>(*original->getGeneBits());
			if (sizeof(float) != sizeof(uint32_t))
				cerr << "float is not " << sizeof(uint32_t) << " bytes long  (same as uint32_t). This could lead to erratic behaviour and loss of information." << endl;
			break;
		case DOUBLE:
			minSeed.doubleValue = original->getMinimumSeed().doubleValue;
			maxSeed.doubleValue = original->getMaximumSeed().doubleValue;
			geneBits = new vector<unsigned char>(*original->getGeneBits());
			if (sizeof(double) != sizeof(uint64_t))
				cerr << "double is not " << sizeof(uint64_t) << " bytes long  (same as uint64_t). This could lead to erratic behaviour and loss of information." << endl;
			break;
		default: //TODO
			break;
	}

}

Gene::~Gene() {
	if (geneBits)
		delete geneBits;
}

inline vector<unsigned char>* Gene::getGeneBits()
{
	return geneBits;
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

inline void Gene::printBits()
{
	if (!geneBits)
		return;
	for (vector<unsigned char>::iterator it = geneBits->begin(); it != geneBits->end(); it++)
		cout << (uint16_t)*it << ", ";
	cout << "\b\b  " << endl;
}

inline void Gene::setValueUInt8(uint8_t value) {	
	int j = 0;
	while (value) {
		(*geneBits)[j] = (value & 1);
		value >>= 1;		
		j++;
	}
}

inline void Gene::setValueInt8(int8_t value) {
	setValueUInt8((uint8_t)value);
}

inline void Gene::setValueUInt16(uint16_t value) {
	int j = 0;
	while (value) {
		(*geneBits)[j] = (value & 1);
		value >>= 1;
		j++;
	}
}

inline void Gene::setValueInt16(int16_t value) {
	setValueUInt16((uint16_t)value);
}

inline void Gene::setValueUInt32(uint32_t value) {
	int j = 0;
	while (value) {
		(*geneBits)[j] = (value & 1);
		value >>= 1;
		j++;
	}
}

inline void Gene::setValueInt32(int32_t value) {
	setValueUInt32((uint32_t)value);
}

inline void Gene::setValueUInt64(uint64_t value) {
	int j = 0;
	while (value) {
		(*geneBits)[j] = (value & 1ui64);
		value >>= 1;
		j++;
	}
}

inline void Gene::setValueInt64(int64_t value) {
	setValueUInt64((uint64_t)value);
}

inline void Gene::setValueFloat(float value) {
	if (sizeof(float) != sizeof(uint32_t))
		cerr << "float is not " << sizeof(uint32_t) << " bytes long  (same as uint32_t). This could lead to erratic behaviour and loss of information." << endl;
	setValueUInt32(*(uint32_t *)(&value));
}

inline void Gene::setValueDouble(double value) {
	if (sizeof(double) != sizeof(uint64_t))
		cerr << ERR_DOUBLE_SIZEMISMATCH << endl;
	setValueUInt64(*(uint64_t *)(&value));
}