// copyright 2016 john howard (orthopteroid@gmail.com)
// MIT license
//
// Genetic Algorithim Maximizer
// Uses sampling from discrete 'breathing' byte-distributions for noise-generation and jump-mutation.
//
// kudos to Andrew Schwartzmeyer for the jump mutation insight

#ifndef PSYCHICSNIFFLE_SNIFFLE_H
#define PSYCHICSNIFFLE_SNIFFLE_H

#include <iostream>
#include <cmath>
#include <cstring>
#include <omp.h>
#include <functional>
#include <assert.h>

#include "taus88.h"

//////////////////////////////////

namespace sniffle {

// performs random selections, without replacement
// uses a small exclude-list to hopefully increase cache-line usage
template<uint M>
struct NSelector
{
    uint N;
    uint n[M];
    uint selected;

    NSelector(uint N) {
        this->N = N;
        selected = 0;
    }

    void reset() { selected = 0; }

    uint select(Taus88 &fnRand)
    {
        if(selected == M)
            throw new std::runtime_error("too many selections");

        uint a = fnRand() % (N - selected);
        for(int i=0; i<selected; i++)
        {
            if(a == n[i]) a = (a + 1) % (N - selected); // skip over holes
        }
        n[selected++] = a;
        return a;
    }
};

// Performs a normalization and a prefix sum on the input array.
// Input and output arrays can use mismatched numeric types.
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

// Splices two byte arrays of the same length together at the specified bit.
// When splicing, bits from a will be written to lower memory than bits from b.
// In the transition byte, high bits come from a and high bits come from b. Does it matter?
template<uint Size>
void splice(uint8_t *out, uint8_t *a, uint8_t *b, uint uRand) {
    const uint8_t u8Mask[] = {
            0b00000000, 0b00000001, 0b00000011, 0b00000111,
            0b00001111, 0b00011111, 0b00111111, 0b01111111
    };

    uint uBit = uRand % ( Size * 8 );
    int i = 0;
    while (i < (uBit / 8)) out[i++] = a[i];
    out[i++] = (a[i] & ~u8Mask[uBit & 7]) | (b[i] & u8Mask[uBit & 7]);
    while (i < Size) out[i++] = b[i];
}

//////////////////////////////////

template<typename StateType>
struct ByteAnalyser {
    const static int StateSize = sizeof(StateType);

    uint8_t *GetByteArr(StateType &state) { return (uint8_t *) &state; }

    uint8_t distr[StateSize][256];
    uint8_t dSampler[StateSize][65535];
    uint16_t dSamplerN[StateSize];

    void dumpStats() {
        for (int ss = 0; ss < StateSize; ss++) {
            for (int b = 0; b < 256; b++) putchar('A' + 25 * (distr[ss][b]) / 255);
            putchar('\n');
        }
        putchar('\n');
    }

    void crank(StateType *stateArr, int *eliteArr, const int eliteSamples) {
#if 0
        dumpStats();
#endif

        // attenuate - BREATHE OUT
#pragma omp parallel for
        for (int ss = 0; ss < StateSize; ss++)
            for (int b = 0; b < 256; b++)
                if (distr[ss][b] > 1) distr[ss][b] -= 1;

        // amplify by sampling from the elite group - BREATHE IN
#pragma omp parallel for
        for (int i = 0; i < eliteSamples; i++) {
            for (int ss = 0; ss < StateSize; ss++) {
                int b = GetByteArr(stateArr[eliteArr[i]])[ss];
                if (distr[ss][b] < 250) distr[ss][b] += 5;
                // slightly distribute locality
                if (b > 0 && distr[ss][b - 1] < 250) distr[ss][b - 1] += 1;
                if (b < 255 && distr[ss][b + 1] < 250) distr[ss][b + 1] += 1;
            }
        }

        // recalc
#pragma omp parallel for
        for (int ss = 0; ss < StateSize; ss++) {
            dSamplerN[ss] = buildSamplerTable<uint8_t, 65535, uint8_t, 256>(&(dSampler[ss][0]), &(distr[ss][0]));
        }

    }

    void reset() {
        // The value used here to initialize the distribution may have an effect
        // on the convergence rate. A higher value will require more attenuation
        // cycles before the S/R ratio get stronger.
        memset(distr, UINT8_MAX >> 2, StateSize * 256); // a uniform distribution

#pragma omp parallel for
        for (int ss = 0; ss < StateSize; ss++) {
            for (int i = 0; i < 256; i++)
                dSampler[ss][i] = i; // a uniform distribution
            dSamplerN[ss] = ~0;
        }
    }

    void mutatebyte(uint8_t *p, Taus88& fnRand) {
        int byte = fnRand() % StateSize;
        p[byte] = dSampler[byte][ fnRand() % dSamplerN[byte] ];
    }

    void randomize(uint8_t *p, Taus88& fnRand) {
        for (int ss = 0; ss < StateSize; ss++) {
            p[ss] = dSampler[ss][ fnRand() % dSamplerN[ss] ];
        }
    }
};

//////////////////////////////////

template<typename StateType, uint Population>
struct Maximizer {
    const static int StateSize = sizeof(StateType);

    const static int Group2End = Population * .30;
    const static int Group3End = Population * .50;
    const static int Group4End = Population * .70;
    const static int Group5End = Population * .80;
    const static int Group6End = Population * .90;

    const static int EliteSamples = 5 + Group3End * .05;
    int eliteSamples[EliteSamples];

    int pa, pb;

    StateType state[2][Population];

    float_t e[Population];
    uint16_t eSampler[65535];
    uint16_t eSamplerN;

    ByteAnalyser<StateType> byteAnalyser;
    Taus88State taus88State;

    uint8_t *GetByteArr(StateType &state) { return (uint8_t *) &state; }

    StateType *GetStateArr() { return &(state[pa][0]); }

    uint8_t *newPop(int i = 0) { return GetByteArr(state[pb][i]); }

    uint8_t *oldPop(int i = 0) { return GetByteArr(state[pa][i]); }

    Maximizer() {
        taus88State.seed();
    }

    void dumpStats() {
        for (int ss = 0; ss < StateSize; ss++) {
            for (int p = 0; p < Population; p++)
                putchar('A' + 25 * oldPop(p)[ss] / 255);
            putchar('\n');
        }
        putchar('\n');
    }

    void reset() {
        pa = 0;
        pb = 1;
        byteAnalyser.reset();

#pragma omp parallel
        {
            Taus88 taus88(taus88State);
#pragma omp for
            for (int i = 0; i < Population; i++) {
                byteAnalyser.randomize(oldPop(i), taus88);
            }
        }
    }

    void crank(float *f) {
#if 0
        dumpStats();
#endif

        // find max, and clobber it to prevent saturation
        int imax = 0;
        int imin = 0;
        for (int i = 1; i < Population; i++) {
            if (f[i] > f[imax]) imax = i;
            if (f[i] < f[imin]) imin = i;
        }
        for (int i = 1; i < Population; i++) {
            if (f[i] == f[imax]) f[i] = f[imin];
        }

        // calc sampler table
        eSamplerN = buildSamplerTable<uint16_t, 65535, float_t, Population>(eSampler, f);

        // sample elites
        eliteSamples[0] = imax; // add best only once to prevent saturation
#pragma omp parallel
        {
            Taus88 taus88(taus88State);
#pragma omp for
            for (int i = 1; i < EliteSamples; i++)
                eliteSamples[i] = eSampler[ taus88() % eSamplerN ];
        }

        // refresh distributions
        byteAnalyser.crank(GetStateArr(), eliteSamples, EliteSamples);

        /////////////////////////////
        // build next generation

        memcpy(newPop(0), oldPop(imax), (uint) StateSize);

#pragma omp parallel
        {
            Taus88 taus88(taus88State);
            NSelector<2> nselector( eSamplerN );

#pragma omp for nowait
            for (int i = 1; i < Group2End; i++) {
                // g2: preserve elites
                int p = eSampler[ taus88() % eSamplerN ];
                memcpy(newPop(i), oldPop(p), (uint) StateSize);
            }
#pragma omp for nowait
            for (int i = Group2End +1; i < Group3End; i++) {
                // g3: semi-preserve elites
                int p = eSampler[ taus88() % eSamplerN ];
                memcpy(newPop(i), oldPop(p), (uint) StateSize);
                byteAnalyser.mutatebyte(newPop(i), taus88);
            }
#pragma omp for nowait
            for (int i = Group3End +1; i < Group4End; i++) {
                // g4: some favourables are spliced with best
                int b = eSampler[ taus88() % eSamplerN ];
                splice<StateSize>(newPop(i), oldPop(0), oldPop(b), (uint) taus88());
            }
#pragma omp for nowait
            for (int i = Group4End +1; i < Group5End; i++) {
                // g5: some favourables are spliced with best (other way)
                int a = eSampler[ taus88() % eSamplerN ];
                splice<StateSize>(newPop(i), oldPop(a), oldPop(0), (uint) taus88());
            }
#pragma omp for nowait
            for (int i = Group5End +1; i < Group6End; i++) {
                nselector.reset();
                // g6: favourables that are only spliced
                int a = eSampler[ nselector.select(taus88) ];
                int b = eSampler[ nselector.select(taus88) ];
                splice<StateSize>(newPop(i), oldPop(a), oldPop(b), (uint) taus88());
            }
#pragma omp for nowait
            for (int i = Group6End +1; i < Population; i++) {
                // g7: randomize rest using byteAnalyser
                byteAnalyser.randomize(newPop(i), taus88);
            }
        }

        std::swap(pa, pb);
    }

};

}

#endif //PSYCHICSNIFFLE_SNIFFLE_H
