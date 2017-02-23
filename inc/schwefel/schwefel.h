#ifndef PSYCHICSNIFFLE_SCHWEFEL_H
#define PSYCHICSNIFFLE_SCHWEFEL_H

//
// Created by orthopteroid on 02/12/16.
//

#include <iostream>
#include <cmath>
#include <cstring>

#include "sniffle.h"

using namespace sniffle;

//////////////////////////////

template<typename Rep, uint Dimension, uint Population>
struct Schwefel
{
    typedef Rep StateType[Dimension];

    static double_t Eval(const StateType& state)
    {
        // schwefel per https://www.sfu.ca/~ssurjano/Code/schwefm.html
        // global min for D=20 is (420.9687, ..., 420.9687) and can be found using best methods in 40K evals
        double_t sum = 0.f;
        for(int d=0; d<Dimension; d++)
        {
            double_t x = (double_t)1000. * (double_t) state[d] / (double_t)((Rep)~0) - (double_t)500.;
            sum += x * sin(sqrt(abs(x)));
        }
        return sum - 418.9829 * Dimension;
    }

    static void Solve(int solns)
    {
        printf("Maximize -Schwefel<%d> : https://www.sfu.ca/~ssurjano/schwef.html\n", Dimension);

        float_t f[Population];
        Maximizer<StateType, Population> solver;

        float_t best = 0.f;
        uint i = 0;

        uint t = 0;
        solver.reset();
        while( t < 1e6 ) {
#pragma omp for
            for( int p=0; p<Population; p++ )
                f[p] = Eval( solver.GetStateArr()[p] );

            if( t == 10000 || best != f[0] || best > 0.f ) {
                best = f[0];
                if( t == 10000 || best > 0. ) { // terminate when +ve or when we hit 10k iterations
                    //solver.byteAnalyser.dumpStats(); // enable this to see internal stats
                    printf("%u, %8.4f, ", t, best);
                    for( int d=0; d<Dimension; d++ ) {
                        StateType &state = *solver.GetStateArr();
                        double_t val = (float_t) 1000. * (float_t) state[d] / (float_t) ((Rep) ~0) - (float_t) 500.;
                        printf("%8.8f%c ", val, d < Dimension - 1 ? ',' : '\n' );
                    }
                    fflush(stdout);
                    t = 0;
                    solver.reset();
                    if( ++i == solns ) break;
                    continue;
                }
            }

            solver.crank(f);
            t++;
        }
    }
};

#endif //PSYCHICSNIFFLE_SCHWEFEL_H
