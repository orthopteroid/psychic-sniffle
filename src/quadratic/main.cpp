// copyright 2018 john howard (orthopteroid@gmail.com)
// MIT license

#include <iostream>
#include <cstring>
#include <time.h>
#include <cstring>
#include <limits>

#include "sniffle.h"

#include "cpuinfo.h"

using namespace sniffle;

//////////////////////////////

timespec diff(timespec start, timespec end)
{
    timespec temp;
    if (end.tv_nsec < start.tv_nsec) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

template<uint Population>
struct Quadratic
{
    float_t f[Population];
    Maximizer<float, Population, ByteAnalyser> solver;

    using Rep = float;

    constexpr static float mu = 101.10101f; // target

    static float_t Eval(const float& x)
    {
        if( isnanf(x) ) return 0.f;
        const float delta = 20.f;
        return (1.f/sqrtf(2.f*delta*delta*(float)M_PI))*expf(-(x-mu)*(x-mu)/(2.f*delta*delta));
    }

    void Solve()
    {
        for( int p=0; p<Population; p++ )
            solver.GetStateArr()[p] = ((float)((rand() % 20000))) - 10000;

        uint iterations = 1;
        solver.reset();

        timespec time1, time2;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);

        while( true ) {
#pragma omp for
            for( int p=0; p<Population; p++ )
                f[p] = Eval( solver.GetStateArr()[p] );

            solver.crank(f);

            float percent = 100 * fabs(mu - *solver.GetStateArr()) / mu;
            //printf("%8.8f, ", percent );
            if( percent < .01 ) break;
            iterations++;
        }

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);

        printf("%8.8f, ", *solver.GetStateArr() );
        if(diff(time1,time2).tv_sec == 0)
            std::cout<<diff(time1,time2).tv_nsec / 1e6 / iterations << ", " << iterations << ", OK\n";
        else
            std::cout << "0, 0, ERR\n";
    }
};

int main()
{
    srand(int(time(NULL)));

#if defined(NDEBUG)
    {
        int cores = EnumCores();
        printf("OpenMP using %d threads\n", cores);
        omp_set_num_threads(cores);
    }
#endif

    printf("Maximize Quadratic. Target = %f\n", Quadratic<60000>::mu);

    for(int i=0;i<10;i++)
    {
        {
            Quadratic<60000> solver;
            solver.Solve();
        }
    }

    return 0;
}
