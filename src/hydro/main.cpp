#include <iostream>
#include <cstring>

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <time.h>
#include <bits/siginfo.h>

#include "hydro/hydro.h"

#include "sniffle.h"

#define ENABLE_PAGINATED_OUTPUT

////////////

using namespace std;
using namespace sniffle;
using namespace hydro;

// hill curve point samples: head vs power vs efficiency
// http://encyclopedia2.thefreedictionary.com/Hydroturbine (fig 6)
float m_pUnitPHE[] =
        {
            20,12, 130,28, /* min p,h  max p,h */
            /* h, p, e, p, e, ... , -1, h, p, e, p, e, ... , -1, h, ... , -1, -1 */
            28, 20,89, 80,93.5, 130,92, -1,
            26, 35,91, 45,92, 50,93, 60,93.5, 80,93.5, 105,93.5, 115,93, 125,92, -1,
            24, 25,90, 35,91, 42,92, 55,93, 60,93.5, 70,93.5, 95,93.5, 105,93, 125,92, -1,
            22, 30,90, 35,91, 42,92, 55,93, 65,93.5, 70,93.5, 90,93, 95,92, 105,91, 120,90, -1,
            21, 65,93.5, -1,
            20, 30,90, 38,91, 45,92, 65,93, 80,92, 90,91, 105,90, 110,89, 130,84, -1,
            18.5, 55,92, -1,
            18, 28,89, 35,90, 42,91, 65,91, 80,90, 85,89, 90,90, 105,84, -1,
            17.5, 50,91, -1,
            16.1, 40,90, -1,
            16, 25,88, 30,89, 45,90, 65,89, 70,88, 85,84, -1,
            15, 35,89, 1,
            14, 35,88, 50,88, 65,84, -1,
            13.5, 35,88, -1,
            12, 20,0, 120,0, -1
                             -1
        };

float m_pFeasZone[] =
        {
            50, /* avg p */
            /* p,h, p,h, ... p,h, -1 */
            41,27, 130,27, 125,22.5, 100,18, 40,13.5, 35,14, 25,17, 35,26, -1 // clockwise
        };

float m_pRoughZone[] =
        {
            /* p,h, p,h, ... p,h, -1 */
            90,15, 60,20, 100,23, -1, // clockwise
        };

float demand[] = { /*0.f,*/ 0.f /* allow no-cost warmup*/, 90.f, 110.f, 90.f, 80.f, 150.f, 210.f, 180.f, 110.f, 90.f, 80.f, 90.f };

const float inQK = 1.5f;
const float inflow[] = { inQK*40.f, inQK*30.f, inQK*20.f, inQK*40.f, inQK*50.f, inQK*40.f, inQK*30.f, inQK*20.f, inQK*40.f, inQK*50.f, inQK*40.f, inQK*30.f };

////////////////

// configuration constants for basin are lumped into a single struct
struct RiverConfig
{
    PlantCoefs<1> upperC; // 1 unit
    PlantCoefs<2> lowerC; // 2 units
};

// river unit-operations are the solver's "decision variables"
// the solver makes guesses about "configurations of river operations over the timescale"
struct RiverOp
{
    UnitOpArr<1> upperU; // 1 unit
    UnitOpArr<2> lowerU; // 2 units
};

template<size_t StepCount>
using RiverOpArr = RiverOp[StepCount];

// river simulations are required to determine the value of guessed river-operations
// the timeseries output for each timestep for the whole basin is lumped together
struct RiverStep
{
    PlantStep<1> upperP;
    UnitStepArr<1> upperU;

    PlantStep<2> lowerP;
    UnitStepArr<2> lowerU;
};

template<size_t StepCount>
using RiverStepArr = RiverStep[StepCount];

// a helper method that initializes the basin's "current state"
void Initialize(RiverStep &s, const RiverConfig& c)
{
    auto fnUnit = [](UnitStep& u) -> void
    {
        u.m_AvgQ = u.m_AvgP = u.m_AvgE = 0;
        u.m_CurState = StateType::STOP;
    };

    // reset upper plant and units
    s.upperP.initialize( c.upperC ); // sadly, no polymophic lambdas in c++11...
    fnUnit( s.upperU[0] );

    // reset lower plant and units
    s.lowerP.initialize( c.lowerC ); // sadly, no polymophic lambdas in c++11...
    fnUnit( s.lowerU[0] );
    fnUnit( s.lowerU[1] );
}

// the simulation & objective function routine for the basin.
// taking an array used to drive the simulation decisions and an array to output the timeseries results.
// also outputs the objective function value to be used by the solver to weigh the simulation's value.
template<uint StepCount>
float Simulate(RiverStepArr<StepCount> &steps, RiverOpArr<StepCount> &ops)
{
    // system integration and mass-conversion coefficients
    const float IMPERIAL = 62.4f /* POUNDSPERCUBICFT */ * 0.746f /* KWPERHP */ / 550.f /* FTPOUNDSPERHP */; /* for cfs from kw */
    const float METRIC   = 1000.0f /* WATERDENSITYINKGPERM3 */ * 9.81f /* ACCELDUETOGRAVITY */ / 1000.f /* WATTSPERKW */; /* for cms from kw */
    SystemCoefs syscoefs =
        {
            1.f / 12.f, // discharge integration coef: storage in xHOURS, discharge in AVGx for 5 min
            IMPERIAL, // power conversion coef
        };

    const RiverConfig conf =
        {
            {
                // upper plant config
                .05f, 12.f, 28.f, .001f, // m_SSslope, m_SSmin, m_SSmax, m_TWslope
                10.f, 20.f, // m_WarmupQ, m_SpinQ
                syscoefs, // reference to shared object
                // upper plant unit
                { m_pUnitPHE },
                { m_pFeasZone },
                { m_pRoughZone },
            },
            {
                // lower plant config
                .025f, 12.f, 28.f, .002f, // m_SSslope, m_SSmin, m_SSmax, m_TWslope
                10.f, 20.f, // m_WarmupQ, m_SpinQ
                syscoefs, // reference to shared object
                // lower plant units
                {m_pUnitPHE,m_pUnitPHE},
                {m_pFeasZone,m_pFeasZone},
                {m_pRoughZone,m_pRoughZone},
            },
        };

    // set initial reservoir storage
    RiverStep initRS;
    Initialize( initRS, conf );

    // simulate
    for( uint t=0; t<StepCount; t++ )
    {
        const RiverStep &prevRS = (t == 0) ? initRS : steps[t-1];

        steps[t].upperP.simulate(
            steps[t].upperU, ops[t].upperU, // current timestep for unit state (output) according to unit operations (input)
            inflow[t],                      // inflow to upper plant comes from data array (input)
            prevRS.upperP, prevRS.upperU,   // previous timestep for upper plant and its units (input)
            conf.upperC
        );

        steps[t].lowerP.simulate(
            steps[t].lowerU, ops[t].lowerU, // current timestep for unit state (output) according to unit operations (input)
            steps[t].upperP.totQS(),        // inflow to lower plant comes from upper plant (input)
            prevRS.lowerP, prevRS.lowerU,   // previous timestep for lower plant and its units (input)
            conf.lowerC
        );
    }

    // collect stats for objective function
    uint ancillary = 0;
    uint roughZone = 0;
    uint starts = 0, stops = 0;
    StatAvg statEff;
    StatPosNeg statPow;
    for( uint t=0; t<StepCount; t++ )
    {
        // upper plant stats
        for( int u = 0; u < conf.upperC.GetUnitCount(); u++ )
        {
            const UnitStep &prev = (t == 0) ? initRS.upperU[u] : steps[t-1].upperU[u];
            const UnitStep &unit = steps[t].upperU[u];
            if( unit.isStarting( prev.m_CurState ) ) starts++;
            if( unit.isStopping( prev.m_CurState ) ) stops++;
            if( unit.isRoughZone( conf.upperC.m_RoughZoneArr[u], steps[t].upperP.m_Head ) ) roughZone++;
        }
        statEff.incGZ( steps[t].upperP.m_AvgE );

        // lower plant stats
        for( int u = 0; u < conf.lowerC.GetUnitCount(); u++ )
        {
            const UnitStep &prev = (t == 0) ? initRS.lowerU[u] : steps[t-1].lowerU[u];
            const UnitStep &unit = steps[t].lowerU[u];
            if( unit.isStarting( prev.m_CurState ) ) starts++;
            if( unit.isStopping( prev.m_CurState ) ) stops++;
            if( unit.isRoughZone( conf.lowerC.m_RoughZoneArr[u], steps[t].lowerP.m_Head ) ) roughZone++;
        }
        statEff.incGZ( steps[t].lowerP.m_AvgE );

        // tally off-demand production
        float pow = steps[t].upperP.m_AvgP + steps[t].lowerP.m_AvgP;
        statPow.inc( pow - demand[t] );
    }
    float powDev = statPow.pos + statPow.neg;
    float totSS = starts + stops;
    float volPosDevUpper = - min( steps[StepCount-1].upperP.m_Vol - initRS.upperP.m_Vol, 0.f ); // no penalty for +ve
    float volPosDevLower = - min( steps[StepCount-1].lowerP.m_Vol - initRS.lowerP.m_Vol, 0.f ); // no penalty for +ve

    // objective function...
    return
//        (1e6f - volPosDevUpper) + (1e6f - volPosDevLower) + // minimize start-to-end pool volume deviation
        (1e6f - powDev) + // minimize deviation from demand
        1e3f * statEff.avg() + // maximize efficiency
        (1e3f - totSS) + // minimize total starts and stops
        (1e3f - roughZone) // minimize roughzone operation
    ;
}

///////////////////////

struct TimerSignaller
{
    timer_t timer;
    TimerSignaller( uint msec, uint signo, uint extra )
    {
        struct sigevent sigev;
        sigev.sigev_notify = SIGEV_SIGNAL;
        sigev.sigev_signo = signo;
        sigev.sigev_value.sival_int = extra;
        if (timer_create(CLOCK_REALTIME, &sigev, &timer) != 0)
        {
            perror("timer_create error!");
        }

        struct itimerspec itval;
        itval.it_value.tv_sec = 1;
        itval.it_value.tv_nsec = 0;//1000000000;// msec * 1000000L;
        itval.it_interval.tv_sec = itval.it_value.tv_sec;
        itval.it_interval.tv_nsec = itval.it_value.tv_nsec;
        if (timer_settime(timer, 0, &itval, 0) != 0)
        {
            perror("time_settime error!");
        }
    }
    ~TimerSignaller()
    {
        struct itimerspec itval;
        itval.it_value.tv_sec = 0;
        itval.it_value.tv_nsec = 0;
        itval.it_interval.tv_sec = itval.it_value.tv_sec;
        itval.it_interval.tv_nsec = itval.it_value.tv_nsec;
        timer_settime(timer, 0, &itval, 0);
        timer_delete(timer);
    }
};

///////////////////////

const uint Steps = 12;
const uint Population = 500;

int lastSignal = 0;
void sig_handler(int sig, siginfo_t * info, void *context)
{
    lastSignal = sig;
    //info->_sifields._timer.si_sigval
}

int main()
{
    srand(int(time(NULL)));

    struct sigaction sigact;
    sigemptyset(&sigact.sa_mask);
    sigact.sa_flags = 0;
    sigact.sa_sigaction = sig_handler;
    sigaction(SIGINT, &sigact, nullptr);

    sigaction(SIGUSR1, &sigact, nullptr);
    //TimerSignaller timerSignaller( 500, SIGUSR1, 42 ); // borked

    printf("multireservoir hydro operations\n");

    // the solver's decision variables are "unit operations" for both reservoirs over the timescale
    Maximizer<RiverOpArr<Steps>, Population> solver;

    // we perform simulations for all the solver's selected unit operations.
    // this is working storage for those simulations.
    RiverStepArr<Steps> steps;

    // each simulation results in an objective value that is then fed back to the solver to tune it's guesses
    float_t f[Population];

    float_t best = -HUGE_VALF; // solver is a maximizer so initialize to -huge_val
    uint iter = 0;
    solver.reset();
    while( true )
    {
        // Simulate the river system using the solver's guesses at what good operations might look like.
        // The solver's convention is that the first guess ( f[0] ) is the "current best guess", so
        // we simulate the guesses is decreasing order so that the last one simulated is the best - which
        // we shall then print.
        for( int p=Population-1; p>=0; p-- ) // note: backwards loop
            f[p] = Simulate<Steps>( steps, solver.GetStateArr()[p] );

        bool terminate = (iter == 10000 || lastSignal == SIGINT );
        if( terminate || f[0] > best || lastSignal ==  SIGUSR1 )
        {
            best = f[0];

            // calc summary stats
            StatAvg statPow, statEff;
            StatMinMax statMMPow;
            for( uint t=0; t<Steps; t++ )
            {
                float stepPow = steps[t].upperP.m_AvgP + steps[t].lowerP.m_AvgP;
                statPow.inc( stepPow );
                statMMPow.inc( stepPow - demand[t] );
                statEff.incGZ( steps[t].upperP.m_AvgE );
                statEff.incGZ( steps[t].lowerP.m_AvgE );
            }

#if defined(ENABLE_PAGINATED_OUTPUT)
            printf("\033c");
#endif

            const char cUnitState[] = {'_','d','g','s','G','S'};

            printf("%5s %5s %6s %2s ", "Qi", "D", "dP", "A" );
            printf("! %6s %5s %5s %5s %5s ", "V", "H", "P", "Q", "S" );
            printf("( %1s %5s %5s ) ", "?", "E", "P" );
            printf("! %6s %5s %5s %5s %5s ", "V", "H", "P", "Q", "S" );
            printf("( %1s %5s %5s %1s %5s %5s ) ", "?", "E", "P", "?", "E", "P" );
            putchar('\n');
            for( uint t=0; t<Steps; t++ )
            {
                printf("%5.1f %5.1f %6.1f %2d ",
                       inflow[t], demand[t], (steps[t].upperP.m_AvgP + steps[t].lowerP.m_AvgP) - demand[t],
                       steps[t].upperP.m_iAncillary + steps[t].lowerP.m_iAncillary
                );
                printf("! %6.1f %5.1f %5.1f %5.1f %5.1f ",
                       steps[t].upperP.m_Vol, steps[t].upperP.m_Head,
                       steps[t].upperP.m_AvgP, steps[t].upperP.m_AvgQ, steps[t].upperP.m_AvgS
                );
                printf("( %c %5.1f %5.1f ) ",
                       cUnitState[ steps[t].upperU[0].getState() ], steps[t].upperU[0].m_AvgE, steps[t].upperU[0].m_AvgP
                );
                printf("! %6.1f %5.1f %5.1f %5.1f %5.1f ",
                       steps[t].lowerP.m_Vol, steps[t].lowerP.m_Head,
                       steps[t].lowerP.m_AvgP, steps[t].lowerP.m_AvgQ, steps[t].lowerP.m_AvgS
                );
                printf("( %c %5.1f %5.1f %c %5.1f %5.1f ) ",
                       cUnitState[ steps[t].lowerU[0].getState() ], steps[t].lowerU[0].m_AvgE, steps[t].lowerU[0].m_AvgP,
                       cUnitState[ steps[t].lowerU[1].getState() ], steps[t].lowerU[1].m_AvgE, steps[t].lowerU[1].m_AvgP
                );
                putchar('\n');
            }
            printf("I %5d E %5.1f P %5.1f MMP %5.1f \n", iter, statEff.avg(), statPow.avg(), statMMPow.maximum() );
            fflush(stdout);

            if( terminate ) break;
        }

        //if( iter > 1 && !(iter % 100) ) solver.byteAnalyser.reset();

        solver.crank(f);
        iter++;

    }
    solver.reset();

    return 0;
}
