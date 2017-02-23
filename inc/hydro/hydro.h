#ifndef PSYCHICSNIFFLE_HYDRO_H
#define PSYCHICSNIFFLE_HYDRO_H

#include <sys/types.h>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "hydro/math.h"

namespace hydro {

//////////////////
// basics of operating units

// the unit running states are the operating modes of the unit
enum class StateType : uint8_t
{
    STOP = 0,		// halted, requires a warmup to be put in use
    SHUTDOWN,		// shutdown takes a 5 min timestep. no power can be generated but some water is used.
    WARMUP_G,		// warmingup to gen is required for a 5 min timestep from stopped state. no power can be generated but some water is used.
    WARMUP_S,		// warmingup to spin is required for a 5 min timestep from stopped state. no power can be generated but some water is used.
    GENERATE,		// currently producing power. does not require warmup to switch to spin state.
    SPIN,			// currently using water but not producing power. does not require warmup to switch to generating state.
    STATES
};
const uint STATES = (uint)StateType::STATES;

// unit operations tell the unit what to do at the start of the current timestep and what state it will be in at the end of the timestep
enum class OpType : uint8_t
{
    CONTINUE = 0,		// continue with currently scheduled state
    SHUTDOWN_STOP = 1,	// begin shutdown process and stopSignal unit at end of timestep
    WARMUP_GEN = 2,		// if stopped, begin warmup process and start generating at end of timestep, otherwise generate now
    WARMUP_SPIN = 3,	// if stopped, begin warmup process and start spinning at end of timestep, otherwise spin now
    OPS,
};
const uint OPS = (uint)OpType::OPS;

inline StateType CalcNextState(StateType prev, uint8_t op)
{
    // nextstate = curstate X op
    const StateType mxState[STATES][OPS] =
        {
            //OpType::CONTINUE      OP_OpType::SHUTDOWN_STOP   OP_WARMUP_GEN           OP_WARMUP_SPIN
            {StateType::STOP,       StateType::STOP,           StateType::WARMUP_G,    StateType::WARMUP_S},    // StateType::STOP
            {StateType::STOP,       StateType::STOP,           StateType::WARMUP_G,    StateType::WARMUP_S},    // StateType::SHUTDOWN
            {StateType::GENERATE,   StateType::STOP,           StateType::GENERATE,    StateType::SPIN},        // StateType::WARMUP_G
            {StateType::SPIN,       StateType::STOP,           StateType::GENERATE,    StateType::SPIN},        // StateType::WARMUP_S
            {StateType::GENERATE,   StateType::SHUTDOWN,       StateType::GENERATE,    StateType::SPIN},        // StateType::GENERATE
            {StateType::SPIN,       StateType::SHUTDOWN,       StateType::GENERATE,    StateType::SPIN},        // StateType::SPIN
        };
    return mxState[ (uint)prev ][ op ];
}

/////////////////
// Search code

struct UnitOp
{
    const static uint8_t OPBITS = 2;
    const static uint8_t FRACBITS = 8 - OPBITS;
    const static uint8_t FRACNORM = (1 << FRACBITS) -1; // for normalization

    union
    {
        uint8_t val;
        struct
        {
            uint8_t op : OPBITS;	// operational change to unit
            uint8_t frac : FRACBITS;	// parameter used for operational change, if necessary
        };
    };

    void set( OpType op_, uint8_t frac_ ) { op = (uint8_t)op_; frac = frac_; }

    float getFrac() const
    {
        const float fracnorm = 1.f / float( FRACNORM );
        return float( frac ) * fracnorm;
    }
};

template<size_t UnitCount>
using UnitOpArr = UnitOp[UnitCount];

////////////////////
// Simulation code

struct SystemCoefs
{
    float m_QIntegrationCoef;
    float m_PConversionCoef;
};

struct UnitStep
{
    // unit state for current timestep
    float m_AvgP;		// power output
    float m_AvgE;		// unit Eff
    float m_AvgQ;		// discharge rate, calc'd from P and E and P->Q discharge assumption
    StateType m_CurState; // current running state

    float getP() const { return m_CurState == StateType::GENERATE ? m_AvgP : 0.f; }
    float getE() const { return m_CurState == StateType::GENERATE ? m_AvgE : 0.f; }
    float getQ() const { return m_AvgQ; }
    float getSpinQ() const { return m_CurState == StateType::SPIN ? m_AvgQ : 0.f; }
    uint getState() const { return (uint)m_CurState; }

    uint getPn() const { return m_CurState == StateType::GENERATE ? 1 : 0; }
    uint getEn() const { return m_CurState == StateType::GENERATE ? 1 : 0; }

    bool isStarting(StateType prev) const { return prev == StateType::STOP && m_CurState != StateType::STOP; }
    bool isStopping(StateType prev) const { return prev != StateType::STOP && m_CurState == StateType::STOP; }
    bool isAncillary() const { return m_CurState == StateType::GENERATE || m_CurState == StateType::SPIN; }
    bool isRoughZone(const float *pRoughZone, float fHead) const { return CalcContains( pRoughZone, getP(), fHead ); }

    void simulate(
        const UnitOp op,
        float fHead,
        float *pUnitPHE, float *pFeasZone,
        float fWarmupQ, float fSpinQ, float fPConvCoef
    )
    {
        switch( m_CurState )
        {
            case StateType::STOP:
            {
                m_AvgQ = m_AvgP = m_AvgE = 0.f;
                break;
            }
            case StateType::WARMUP_G:
            case StateType::WARMUP_S:
            case StateType::SHUTDOWN:
            {
                // no P or E. Q is warmup Q
                m_AvgP = m_AvgE = 0;
                m_AvgQ = fWarmupQ;
                break;
            }
            case StateType::SPIN:
            {
                // no P or E. Q is spin Q. A is ancillary benefit.
                m_AvgP = m_AvgE = 0.f;
                m_AvgQ = fSpinQ;
                break;
            }
            case StateType::GENERATE:
            {
                // P is calculated form the frac of the op
                float pmin, pspan;
                CalcSpan( pmin, pspan, pFeasZone, fHead );
                if( pmin * pspan < 1.f )
                {
                    throw new runtime_error("damn. can't think of a meaningful thing to say.");
                }
                m_AvgP = pmin + pspan * op.getFrac(); // using new span and new P
                m_AvgE = CalcInterpolate( pUnitPHE, m_AvgP, fHead );
                m_AvgQ = CalcQ( m_AvgP, m_AvgE, fHead, fPConvCoef );
                break;
            }
        }
    }
};

using FloatPtrArr = float*;

template<uint UnitCount>
struct PlantCoefs
{
    float m_SSslope;        // slope of the stage-storage curve
    float m_SSmin, m_SSmax; // min & max of stage storage curve, in linear distance
    float m_TWslope;        // slope fo the tailwater curve

    // same for all units in plant
    float m_WarmupQ;        // warmup or shutdown Q
    float m_SpinQ;          // spin Q. likely larger than warmup Q.

    SystemCoefs &m_SysCoefs; // reference to shared system coef struct

    const uint GetUnitCount() const { return UnitCount; }

    // for each unit at this plant...
    FloatPtrArr m_PHEArr[UnitCount];        // power X head X efficiency surface, stored as a special point array
    FloatPtrArr m_FeasZoneArr[UnitCount];   // polygon describing the feasible region, with leading powAvg value
    FloatPtrArr m_RoughZoneArr[UnitCount];  // polygon describing the roughzone region
};

template<uint UnitCount>
using UnitStepArr = UnitStep[UnitCount];

template<uint UnitCount>
struct PlantStep
{
    // plant state for current timestep
    float m_Vol;  // reservoir volume
    float m_Head; // net-head
    float m_AvgP; // plant power, only for running units
    float m_AvgE; // plant efficiency, only for running units
    float m_AvgQ; // plant power discharge
    float m_AvgS; // reservoir spill, sometimes required for mass continuity
    uint8_t m_iAncillary;

    // the continuty adjustor facilitates an iterative continuity calulation
    struct ContinuityAdjustor
    {
        const PlantCoefs<UnitCount>& coefs;
        const float avgI;

        PlantStep& plantstep;
        StatAvg statQ, statP, statE;
        float adjV, newPondElev;

        ContinuityAdjustor(
            PlantStep& _plantstep, const float _avgI, const PlantCoefs<UnitCount>& _coefs
        )
            : plantstep(_plantstep), coefs(_coefs),
              avgI(_avgI), adjV(0), newPondElev(0)
        {}

        void calc( UnitStepArr<UnitCount>& unitArr, const UnitOpArr<UnitCount>& unitOpArr )
        {
            statQ.clear();
            statP.clear();
            statE.clear();
            for( size_t u = 0; u < UnitCount; u++ )
            {
                unitArr[u].simulate(
                    unitOpArr[u],
                    plantstep.m_Head,
                    coefs.m_PHEArr[u], coefs.m_FeasZoneArr[u],
                    coefs.m_WarmupQ, coefs.m_SpinQ,
                    coefs.m_SysCoefs.m_PConversionCoef
                );
                statQ.incGZ( unitArr[ u ].getQ() );
                statP.incGZ( unitArr[ u ].getP() );
                statE.incGZ( unitArr[ u ].getE() );
            }
            plantstep.m_AvgQ = statQ.tot; // was avg
            plantstep.m_AvgS = 0;

            // calc pond elev and spill, for sake of continuity permit -ve volumes
            adjV = avgI;
            newPondElev = coefs.m_SSslope * ( plantstep.m_Vol + coefs.m_SysCoefs.m_QIntegrationCoef * ( adjV - plantstep.m_AvgQ ) );
            if( newPondElev > coefs.m_SSmax )
            {
                plantstep.m_AvgS = ( newPondElev - coefs.m_SSmax ) / ( coefs.m_SSslope * coefs.m_SysCoefs.m_QIntegrationCoef );
                adjV -= plantstep.m_AvgS;
                newPondElev = coefs.m_SSmax;
            }
            else if( newPondElev < coefs.m_SSmin )
            {
                newPondElev = coefs.m_SSmin;
            }
        }
    };

    float totQS() const { return m_AvgQ + m_AvgS; }

    void initialize(const PlantCoefs<UnitCount>& coefs)
    {
        // reset to half-full
        m_Head = (coefs.m_SSmin + coefs.m_SSmax) / 2.f;
        m_Vol = m_Head / coefs.m_SSslope;
    };

    // no losses: hydraulic, yard, generator
    // common plant pool, common unit tailwater
    void simulate(
        UnitStepArr<UnitCount> &unitArr, const UnitOpArr<UnitCount> &unitOpArr,
        const float avgI,
        const PlantStep prevPlant, const UnitStepArr<UnitCount> &prevUnitArr,
        const PlantCoefs<UnitCount> &coefs
    )
    {
        m_Head = prevPlant.m_Head;
        m_Vol = prevPlant.m_Vol;

        // apply the unit operation to the unit state, for the current timestep
        // unit 'operations' are "deltas" to the current unit state
        for( size_t u = 0; u < UnitCount; u++ )
        {
            unitArr[u].m_CurState = CalcNextState(prevUnitArr[u].m_CurState, unitOpArr[u].op);
        }

        // iterate on unit operation and average the operating head.
        ContinuityAdjustor adjustor( *this, avgI, coefs );

        adjustor.calc( unitArr, unitOpArr );
        m_Head = ( m_Head + adjustor.newPondElev - coefs.m_TWslope * ( m_AvgQ + m_AvgS ) ) / 2.f;

        adjustor.calc( unitArr, unitOpArr );
        m_Head = ( m_Head + adjustor.newPondElev - coefs.m_TWslope * ( m_AvgQ + m_AvgS ) ) / 2.f;

        // update plant stats
        m_Vol = prevPlant.m_Vol + coefs.m_SysCoefs.m_QIntegrationCoef * ( adjustor.adjV - m_AvgQ - m_AvgS );
        m_AvgP = adjustor.statP.tot;
        m_AvgE = adjustor.statE.avg();

        m_iAncillary = 0;
        for( size_t u = 0; u < UnitCount; u++ )
        {
            if( unitArr[u].isAncillary() ) m_iAncillary++;
        }
    }
};

////////////////

void Test();

}

#endif //PSYCHICSNIFFLE_HYDRO_H
