// copyright 2016 john howard (orthopteroid@gmail.com)
// MIT license
//
// Performs a normalization and a prefix sum on the input array.
// Input and output arrays can use mismatched numeric types.

#ifndef PROJECT_SAMPLERTABLE_H
#define PROJECT_SAMPLERTABLE_H

namespace util {

template<typename OT, uint ON, typename IT, uint IN>
int buildSamplerTable(OT* outArr, IT* inArr)
{
    IT min = inArr[0];
    for(int i=1; i<IN; i++) {
        min = std::min(min, inArr[i]);
    }

    // compute area
    float_t sum = 0;
    for(int i=0; i<IN; i++) {
        sum += (inArr[i] - min);
    }

    // when all equal, return a uniform distr
    if (sum == 0) {
        for(int i=0; i<IN; i++)
            outArr[i] = i;
        return IN;
    }

    // select a scaling coef such that output outArr is filled
    float_t coef = (float)(ON -1) / sum;

    // otherwise return a nonuniform distribution
    uint i = 0;
    for (uint j = 0; j < IN; j++) {
        uint expectedValue = ( (float_t)(inArr[j]) - (float_t)min ) * coef;
        for (uint k = 0; k < expectedValue; k++)
            outArr[i++] = j;
    }

    if (i == 0 || i > ON)
        throw new std::runtime_error("outArr error");

    return i;
}

}

#endif //PROJECT_SAMPLERTABLE_H
