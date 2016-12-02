
#ifndef PSYCHICSNIFFLE_SNIFFLE_H
#define PSYCHICSNIFFLE_SNIFFLE_H

//
// Genetic Algorithim Maximizer
// Uses sampling from discrete 'breathing' byte-distributions for noise-generation and jump-mutation.
//
// Created by orthopteroid on 30/11/16.
// kudos to Andrew Schwartzmeyer

#include <iostream>
#include <cmath>
#include <cstring>
#include <omp.h>
#include <functional>

#include "taus88.h"

//////////////////////////////////

namespace sniffle {

template<typename T>
uint sample(T *arr, int samples, uint32_t randValue) {
    T target = arr[samples - 1] * ((float_t) randValue / (float_t) UINT32_MAX);
    for (int i = 0; i < samples; i++)
        if (arr[i] >= target) return i;
    throw new std::runtime_error("sample overflow");
    return 0;
}

template<typename T1, typename T2>
void accumulate(T1 *arrOut, T2 *arr, int samples) {
    T1 sum = 0;
    for (int i = 0; i < samples; i++)
        arrOut[i] = (sum += arr[i]);
}

void splice(uint8_t *out, uint8_t *a, uint8_t *b, uint uSize, uint uBit) {
    const uint8_t u8Mask[] = {
            0b00000000, 0b00000001, 0b00000011, 0b00000111,
            0b00001111, 0b00011111, 0b00111111, 0b01111111
    };

    int i = 0;
    while (i < (uBit / 8)) out[i++] = a[i];
    out[i++] = (a[i] & ~u8Mask[uBit & 7]) | (b[i] & u8Mask[uBit & 7]);
    while (i < uSize) out[i++] = b[i];
}

//////////////////////////////////

template<typename StateType>
struct ByteAnalyser {
    const static int StateSize = sizeof(StateType);

    uint8_t *GetByteArr(StateType &state) { return (uint8_t *) &state; }

    uint8_t distr[StateSize][256];
    uint32_t dSum[StateSize][256];

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
#pragma omp for
        for (int ss = 0; ss < StateSize; ss++)
            for (int b = 0; b < 256; b++)
                if (distr[ss][b] > 1) distr[ss][b] -= 1;

        // amplify by sampling from the elite group - BREATHE IN
#pragma omp for
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
#pragma omp for
        for (int ss = 0; ss < StateSize; ss++)
            accumulate(&(dSum[ss][0]), &(distr[ss][0]), 256);
    }

    void reset() {
        // The value used here to initialize the distribution may have an effect
        // on the convergence rate. A higher value will require more attenuation
        // cycles before the S/R ratio get stronger.
        memset(distr, UINT8_MAX >> 2, StateSize * 256); // a uniform distribution
    }

    void mutatebyte(uint8_t *p, int byte, uint32_t randValue) {
        p[byte] = sample(&(dSum[byte][0]), 256, randValue); // aka jumping mutation
    }

    void randomize(uint8_t *p, std::function<uint32_t()> fnRand) {
        for (int ss = 0; ss < StateSize; ss++)
            p[ss] = sample(&(dSum[ss][0]), 256, fnRand()); // synthesis using the byte-distributions
    }
};

//////////////////////////////////

template<typename StateType, uint Population>
struct Maximizer {
    const static int StateSize = sizeof(StateType);

    const static int Offset_Elite = Population / 4;

    const static int EliteSamples = 5 + Offset_Elite * StateSize / 200; // /200 is .5% mark
    int eliteSamples[EliteSamples];

    int pa, pb;

    StateType state[2][Population];

    float_t fSum[Population];
    float_t e[Population];
    float_t eSum0[Population]; // includes all items, even item 0
    float_t eSum1[Population - 1]; // does not include item 0. eSum1[0] == e[1]

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
#pragma omp parallel for
            for (int i = 0; i < Population; i++) {
                std::function<uint32_t()> fRand = [&taus88]() -> uint32_t { return taus88.operator()(); };

                byteAnalyser.randomize(oldPop(i), fRand);
            }
        }
    }

    void crank(float *f) {
#if 0
        dumpStats();
#endif

        accumulate(fSum, f, Population);

        // calc expected value
#pragma omp for
        for (int p = 0; p < Population; p++)
            e[p] = f[p] / fSum[Population - 1];

        // calc cumulative prob density
        accumulate(eSum0, &(e[0]), Population);
        accumulate(eSum1, &(e[1]), Population - 1); // note: eSum1 doesn't account for item 0

        // sample elites
        eliteSamples[0] = 0; // add best only once to prevent saturation
#pragma omp parallel
        {
            Taus88 taus88(taus88State);
#pragma omp parallel for
            for (int i = 1; i < EliteSamples; i++)
                eliteSamples[i] =
                        sample(eSum1, Offset_Elite - 1, taus88()) + 1; // note: eSum1 doesn't account for item 0
        }

        // refresh distributions
        byteAnalyser.crank(GetStateArr(), eliteSamples, EliteSamples);

        /////////////////////////////
        // build next generation

        // find and keep best
        int imax = 0;
        for (int i = 1; i < Population; i++)
            if (f[i] > f[imax]) imax = i;
        memcpy(newPop(0), oldPop(imax), (uint) StateSize);

#pragma omp parallel
        {
            Taus88 taus88(taus88State);
#pragma omp parallel for
            for (int i = 1; i < Population; i++) {
                std::function<uint32_t()> fRand = [&taus88]() -> uint32_t { return taus88.operator()(); };

                if (i < Offset_Elite * 0.95) // TODO: fix as Population * frac
                {
                    // preserve elites
                    int p = sample(eSum1, Population - 1, taus88()) + 1; // +1 to skip item 0 and avoid saturation
                    memcpy(newPop(i), oldPop(p), (uint) StateSize);
                } else if (i < Offset_Elite) // TODO: fix as Population * frac
                {
                    // semi-preserve elites
                    int p = sample(eSum1, Population - 1, taus88()) + 1; // +1 to skip item 0 and avoid saturation
                    memcpy(newPop(i), oldPop(p), (uint) StateSize);
                    byteAnalyser.mutatebyte(newPop(i), taus88() % StateSize, taus88());
                } else if (i < Population * 0.50) {
                    // some favourables are spliced with best
                    int a = 0;
                    int b = sample(eSum1, Population - 1, taus88()) + 1; // +1 to skip item 0 and avoid saturation
                    splice(newPop(i), oldPop(a), oldPop(b), (uint) StateSize, (uint) taus88() % (8 * StateSize));
                } else if (i < Population * 0.75) {
                    // some favourables are spliced with best
                    int a = sample(eSum1, Population - 1, taus88()) + 1; // +1 to skip item 0 and avoid saturation
                    int b = 0;
                    splice(newPop(i), oldPop(a), oldPop(b), (uint) StateSize, (uint) taus88() % (8 * StateSize));
                } else if (i < Population * 0.95) {
                    // favourables that are only spliced
                    int a = sample(eSum1, Population - 1, taus88()) + 1; // +1 to skip item 0 and avoid saturation
                    int b = sample(eSum1, Population - 1, taus88()) + 1; // +1 to skip item 0 and avoid saturation
                    splice(newPop(i), oldPop(a), oldPop(b), (uint) StateSize, (uint) taus88() % (8 * StateSize));
                } else {
                    // randomize rest using byteAnalyser
                    byteAnalyser.randomize(newPop(i), fRand);
                }
            }
        }

        std::swap(pa, pb);
    }

};

}

#endif //PSYCHICSNIFFLE_SNIFFLE_H
