#include <list>
#include <vector>
#include <bitset>
#include <random>
#include <ctime>

//mt19937 gen(static_cast<unsigned int>(std::time(0)));
//uniform_int_distribution<int32_t> disLong(LONG_MIN, LONG_MAX);
//uniform_int_distribution<int32_t> disBin(0, 1);

//#define random() disLong(gen)
#define randomBinary() 0//disBin(gen)


void crossOver(vector<int32_t> &parentX, vector<int32_t> &parentY){
	/*
	Uniform crossover
	More info at: https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)
	*/
    int32_t vectorSize = parentX.size();
    vector<int32_t> *childX = new vector<int32_t>(vectorSize);
    vector<int32_t> *childY = new vector<int32_t>(vectorSize);

    bool first = false;

    for(int32_t i = 0; i<vectorSize; i++){
        bitset<32> bitX(parentX[i]);
        bitset<32> bitY(parentY[i]);

        for(int j=0; j<32; j++){
            /*
            The first bit is always kept untouched
            This is a premise of crossover
            */
            if ((i==0) && (j==0))
                continue;

            if(randomBinary()){
                /*
                Set bitX[j] as bitY[j] and vice-versa
                */
                bool bit = bitX.test(j);
                bitX.set(j, bitY.test(j));
                bitY.set(j, bit);
            }

        }

        int32_t x = (int32_t)bitX.to_ulong();
        int32_t y = (int32_t)bitY.to_ulong();

        (*childX)[i]=x;
        (*childY)[i]=y;

    }

}
